# Unit-test fixtures for sc-droplets' DropletUtils primitives, from R and DropletUtils itself.
#
#     Rscript oracles/dropletutils/fixtures.R tests/oracles/dropletutils/fixtures
#
# Run with DropletUtils 1.30 (R 4.5). Reads the Monte Carlo probe grids already in the folder
# (mc_grid_*.tsv.gz: total, observed value) and writes, beside them:
#   mc_counts_<alpha>.txt   DropletUtils:::montecarlo_pval counts, 200 iterations, set.seed(100)
#   spline_<n>.tsv.gz          x, y, then smooth.spline(x, y, df = 20)'s fitted values and derivatives
#   spline_<n>_spar.txt     its spar and df
#   goodturing_<k>.tsv.gz      counts, edgeR::goodTuringProportions
#   optimize.txt, quantile.txt, nknots.txt, padjust.txt, lgamma.txt, sum.txt (+ sum_values.txt.gz)
#   platform.txt            R's platform: lgamma and Good-Turing carry its C library's log and exp
suppressMessages(library(DropletUtils))
out <- commandArgs(trailingOnly = TRUE)[1]
fmt <- function(x) sprintf("%.17g", x)
writeLines(R.version$platform, file.path(out, "platform.txt"))  # whose C library the values carry

ambient <- scan(file.path(out, "mc_ambient.txt"), quiet = TRUE)
set.seed(100)
seeds <- dqrng::generateSeedVectors(200, nwords = 2)
for (a in c("inf", "0.5", "50", "2000")) {
  g <- read.delim(gzfile(file.path(out, sprintf("mc_grid_%s.tsv.gz", a))), header = FALSE)
  r <- rle(g[[1]])
  alpha <- if (a == "inf") Inf else as.numeric(a)
  counts <- DropletUtils:::montecarlo_pval(as.integer(r$values), as.integer(r$lengths), g[[2]], ambient, 200L, alpha, seeds, 1:200)
  writeLines(as.character(counts), file.path(out, sprintf("mc_counts_%s.txt", a)))
}

set.seed(7)
for (n in c(30, 120, 700, 2500)) {
  x <- sort(log10(cumsum(sample(1:3, n, replace = TRUE))))
  y <- 4 - 1.5 * plogis((x - median(x)) * 6) + rnorm(n, sd = 0.02)
  fit <- smooth.spline(x, y, df = 20)
  d1 <- predict(fit, deriv = 1)$y; d2 <- predict(fit, deriv = 2)$y
  writeLines(sprintf("%s\t%s\t%s\t%s\t%s", fmt(x), fmt(y), fmt(fit$y), fmt(d1), fmt(d2)), gzfile(file.path(out, sprintf("spline_%d.tsv.gz", n))))
  writeLines(c(fmt(fit$spar), fmt(fit$df)), file.path(out, sprintf("spline_%d_spar.txt", n)))
}

set.seed(11)
cases <- list(rnbinom(3000, mu = 2, size = 0.3), rpois(500, 4), c(rep(0, 50), rpois(400, 30)), rnbinom(20000, mu = 0.5, size = 0.1))
for (k in seq_along(cases)) {
  x <- cases[[k]]
  writeLines(sprintf("%d\t%s", as.integer(x), fmt(edgeR::goodTuringProportions(x))), gzfile(file.path(out, sprintf("goodturing_%d.tsv.gz", k))))
}

f <- function(a) -((a - 53.2) * (a - 53.2)) / (1 + a / 1000)  # basic arithmetic only: no libm in the path
writeLines(c(fmt(optimize(f, c(0.01, 10000), maximum = TRUE)$maximum), fmt(optimize(function(a) (a - 2.5)^2, c(0, 10))$minimum)),
           file.path(out, "optimize.txt"))
q <- c((1:3000 * 7919) %% 10007, rep(5000, 40))
writeLines(c(fmt(quantile(q, 0.99, names = FALSE)), fmt(quantile(rep(c(1, 7), 50), 0.5, names = FALSE)),
             fmt(quantile((1:17)^2 / 3, 0.37, names = FALSE)), fmt(quantile(c(rep(3, 2969), rep(9, 31)), 0.99, names = FALSE))),
           file.path(out, "quantile.txt"))
writeLines(as.character(sapply(c(10, 49, 50, 199, 200, 799, 800, 3199, 3200, 10000), .nknots.smspl)), file.path(out, "nknots.txt"))
lg <- c(1e-310, 0.003, 0.5, 0.9999, 1, 1.5, 2, 2.5, 7.25, 9.999, 10, 10.5, 308.91341196583244, 409.9, 4934721, 1e18,
        (1:2000) / 7.3, 0.001 * 1.37^(0:60))
writeLines(sprintf("%s\t%s", fmt(lg), fmt(lgamma(lg))), file.path(out, "lgamma.txt"))
set.seed(3)
v <- c(1e16, rep(1, 10), -1e16, rnorm(50000, 1460, 30), rnorm(10000, 0, 1e-3))
writeLines(c(fmt(sum(v)), fmt(sum(v[1:12])), fmt(tail(cumsum(v[1:12]), 1))), file.path(out, "sum.txt"))
writeLines(fmt(v), gzfile(file.path(out, "sum_values.txt.gz")))
p <- c(0.01, NA, 0.04, 0.04, 0.2, 0, NA, 0.9, 0.03, 1)
writeLines(fmt(p.adjust(p, "BH")), file.path(out, "padjust.txt"))
cat("fixtures ->", out, "\n")
