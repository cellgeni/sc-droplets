# Run zUMIs' own cellBC (automatic detection) on a dataset's counts, with a given inflection.
#   Rscript oracle.R BARCODEIDFUN_R INFLECTION_DIR COUNTS.tsv OUT_DIR
# BARCODEIDFUN_R is zUMIs' barcodeIDFUN.R; INFLECTION_DIR a checkout of the R package inflection,
# whose R/*.R are sourced (so its version is the one checked out, not an installed one).
# COUNTS.tsv is `barcode<TAB>count` (oracles/make_counts.py). Writes OUT_DIR/cells.tsv
# (`nReadsperCell<TAB>cells<TAB>smallest count`), OUT_DIR/cells_<nReadsperCell>.txt.gz (the
# barcodes, sorted) and OUT_DIR/versions.txt.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(data.table))
for (f in sort(list.files(file.path(args[2], "R"), pattern = "[.]R$", full.names = TRUE))) source(f)
source(args[1])
dir.create(args[4], recursive = TRUE, showWarnings = FALSE)
lines <- character(0)
for (nr in c(100, 10)) {
  opt <<- list(barcodes = list(nReadsperCell = nr))
  cells <- suppressWarnings(cellBC(bcfile = NULL, bcnum = NULL, bcauto = FALSE, bccount_file = args[3], outfilename = NULL))
  lines <- c(lines, sprintf("%d\t%d\t%d", nr, nrow(cells), as.integer(min(cells$n))))
  out <- gzfile(file.path(args[4], sprintf("cells_%d.txt.gz", nr)), "w")
  writeLines(sort(cells$XC, method = "radix"), out)
  close(out)
}
writeLines(lines, file.path(args[4], "cells.tsv"))
writeLines(c(R.version.string, paste("platform", R.version$platform), paste("data.table", packageVersion("data.table"))),
           file.path(args[4], "versions.txt"))
