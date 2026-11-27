# Fit mclust's two-component mixture of unequal variances to log10 UMIs of the barcodes with at
# least 100, to compare with `gmm`'s converged fit.
#   Rscript mclust.R FEATURES.tsv OUT_DIR
# FEATURES.tsv is `barcode<TAB>umis<TAB>genes` (make_features.py). Writes OUT_DIR/mclust.tsv
# (`components<TAB>model<TAB>log-likelihood<TAB>means<TAB>cells`, cells above 0.5 for the component
# with the highest mean) and OUT_DIR/versions.txt. mclust initialises on a random subset of large
# inputs: set.seed(0) makes it repeatable.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(mclust))
features <- read.delim(args[1], header = FALSE, col.names = c("barcode", "umis", "genes"),
                       colClasses = c("character", "numeric", "numeric"))
x <- log10(features$umis[features$umis >= 100])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
set.seed(0)
fit <- Mclust(x, G = 2, modelNames = "V", verbose = FALSE)
top <- which.max(fit$parameters$mean)
line <- sprintf("2\tV\t%.17g\t%s\t%d", fit$loglik, paste(sprintf("%.17g", fit$parameters$mean), collapse = ","),
                sum(fit$z[, top] > 0.5))
cat(line, "\n")
writeLines(line, file.path(args[2], "mclust.tsv"))
writeLines(c(R.version.string, paste("platform", R.version$platform),
             paste("mclust", as.character(packageVersion("mclust")))),
           file.path(args[2], "versions.txt"))
