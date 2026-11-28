# Run dropbead's own estimateCellNumber on a dataset's counts.
#   Rscript oracle.R TOOLS_R COUNTS.tsv OUT_DIR
# TOOLS_R is dropbead's R/Tools.R; COUNTS.tsv is `barcode<TAB>count` (oracles/make_counts.py).
# Writes OUT_DIR/cells.tsv (`max.cells<TAB>cells`) and OUT_DIR/versions.txt.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(methods))
source(args[1])
counts <- read.delim(args[2], header = FALSE, col.names = c("barcode", "n"), colClasses = c("character", "numeric"))
reads <- sort(counts$n, decreasing = TRUE)
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
lines <- character(0)
for (m in c(15000, 10000)) {
  k <- estimateCellNumber(reads, max.cells = m)
  lines <- c(lines, sprintf("%d\t%d", m, k))
  cat(sprintf("max.cells %d: %d cells\n", m, k))
}
writeLines(lines, file.path(args[3], "cells.tsv"))
writeLines(c(R.version.string, paste("platform", R.version$platform)), file.path(args[3], "versions.txt"))
