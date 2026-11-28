# Call scuttle's isOutlier on a dataset's counts above each floor, as the `mad` rule does.
#   Rscript mad.R COUNTS.tsv OUT_DIR
# COUNTS.tsv is `barcode<TAB>count` (oracles/make_counts.py). Writes OUT_DIR/mad.tsv
# (`floor<TAB>nmads<TAB>lower limit<TAB>cells`, the limit as 17 significant digits) and
# OUT_DIR/versions.txt.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(scuttle))
counts <- read.delim(args[1], header = FALSE, col.names = c("barcode", "n"), colClasses = c("character", "numeric"))
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
lines <- character(0)
for (floor in c(50, 100, 200)) {
  for (nmads in c(2, 3)) {
    x <- counts$n[counts$n >= floor]
    outlier <- isOutlier(x, nmads = nmads, type = "lower", log = TRUE)
    lower <- attr(outlier, "thresholds")[["lower"]]
    lines <- c(lines, sprintf("%d\t%d\t%.17g\t%d", floor, nmads, lower, sum(!outlier)))
    cat(sprintf("floor %d, nmads %d: limit %.6g, %d cells\n", floor, nmads, lower, sum(!outlier)))
  }
}
writeLines(lines, file.path(args[2], "mad.tsv"))
writeLines(c(R.version.string, paste("platform", R.version$platform),
             paste("scuttle", as.character(packageVersion("scuttle")))),
           file.path(args[2], "versions.txt"))
