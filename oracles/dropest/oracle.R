# Run dropestr's own EstimateCellsNumber on a dataset's counts.
#   Rscript oracle.R CELL_NUMBER_R COUNTS.tsv OUT_DIR
# CELL_NUMBER_R is dropestr's R/cell_number.R; COUNTS.tsv is `barcode<TAB>count`
# (oracles/make_counts.py). Writes OUT_DIR/bounds.tsv (`lag<TAB>expected<TAB>min<TAB>max`) and
# OUT_DIR/versions.txt.
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
counts <- read.delim(args[2], header = FALSE, col.names = c("barcode", "n"), colClasses = c("character", "numeric"))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
lines <- character(0)
for (lag in c(0.05, 0.1)) {
  r <- EstimateCellsNumber(counts$n, lag = lag)
  lines <- c(lines, sprintf("%s\t%d\t%d\t%d", format(lag), as.integer(r$expected), as.integer(r$min), as.integer(r$max)))
  cat(sprintf("lag %s: expected %d, min %d, max %d\n", format(lag), as.integer(r$expected), as.integer(r$min), as.integer(r$max)))
}
writeLines(lines, file.path(args[3], "bounds.tsv"))
writeLines(c(R.version.string, paste("platform", R.version$platform)), file.path(args[3], "versions.txt"))
