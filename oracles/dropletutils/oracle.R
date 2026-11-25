# Run DropletUtils' cell calling on one raw matrix and write what sc-droplets compares.
#
#     Rscript oracle.R MTX_DIR OUT_DIR [variant ...]
#
# MTX_DIR holds matrix.mtx (features x barcodes), barcodes.tsv and features.tsv
# (oracles/starsolo/make_mtx.py writes them). Only "Gene Expression" rows are used,
# as sc-droplets does by default. Writes, for the installed DropletUtils:
#   barcoderanks.txt          knee and inflection, %.17g
#   defaultdrops.txt          0-based indices of the cells
#   emptydrops.tsv            per tested barcode: index, total, n_above, FDR; header: alpha, retain
#   emptydropscellranger.tsv  the same, DropletUtils 1.14 and later
# Every emptyDrops call follows set.seed(100). Extra arguments name variants that need
# something other than the installed version's defaults:
#   br128     barcodeRanks from barcodeRanks-1.28.R (evaluated inside DropletUtils 1.30),
#             and emptyDrops with that knee as retain
#   alphainf  emptyDrops(alpha = Inf) (1.2, whose default Dirichlet draws use R's rgamma)
suppressMessages({ library(Matrix); library(DropletUtils) })
args <- commandArgs(trailingOnly = TRUE)
mtx <- args[1]; out <- args[2]; extra <- args[-(1:2)]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

m <- as(readMM(file.path(mtx, "matrix.mtx")), "CsparseMatrix")
features <- read.delim(file.path(mtx, "features.tsv"), header = FALSE, colClasses = "character")
m <- m[features[[3]] == "Gene Expression", , drop = FALSE]
version <- as.character(packageVersion("DropletUtils"))
fmt <- function(x) sprintf("%.17g", x)

points <- function(br) {
  if (!is.null(br$knee)) c(br$knee, br$inflection)  # 1.0 - 1.2 return a list
  else c(S4Vectors::metadata(br)$knee, S4Vectors::metadata(br)$inflection)
}

ranks <- if ("br128" %in% extra) {
  env <- new.env(parent = asNamespace("DropletUtils"))
  for (e in parse(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "barcodeRanks-1.28.R")))
    if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
  points(env$.barcode_ranks(m))
} else points(barcodeRanks(m))
writeLines(c(paste("DropletUtils", version, paste(extra, collapse = " ")), paste("knee", fmt(ranks[1])), paste("inflection", fmt(ranks[2]))),
           file.path(out, "barcoderanks.txt"))

writeLines(as.character(which(defaultDrops(m)) - 1L), file.path(out, "defaultdrops.txt"))

write_test <- function(e, file) {
  tested <- which(!is.na(e$PValue))
  meta <- S4Vectors::metadata(e)
  header <- c(paste("# DropletUtils", version),
              paste("# alpha", if (is.null(meta$alpha)) "NA" else fmt(meta$alpha)),
              paste("# retain", if (is.null(meta$retain)) "NA" else fmt(meta$retain)))
  niters <- 10000
  body <- sprintf("%d\t%d\t%d\t%s", tested - 1L, as.integer(e$Total[tested]),
                  as.integer(round(e$PValue[tested] * (niters + 1) - 1)), fmt(e$FDR[tested]))
  writeLines(c(header, body), file)
}

set.seed(100)
e <- if ("br128" %in% extra) emptyDrops(m, retain = ranks[1]) else if ("alphainf" %in% extra) emptyDrops(m, alpha = Inf) else emptyDrops(m)
if ("br128" %in% extra) S4Vectors::metadata(e)$retain <- ranks[1]
write_test(e, file.path(out, "emptydrops.tsv"))

if (exists("emptyDropsCellRanger")) {
  set.seed(100)
  e <- tryCatch(emptyDropsCellRanger(m), error = function(err) conditionMessage(err))
  if (is.character(e)) writeLines(paste("# error", e), file.path(out, "emptydropscellranger.tsv"))
  else write_test(e, file.path(out, "emptydropscellranger.tsv"))
}
cat("DropletUtils", version, "->", out, "\n")
