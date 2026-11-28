# Run DropletQC's own identify_empty_drops (sourced from the tag in run.sh) on every input table, and
# write what it decided: the cut and ks's bandwidth (summary.tsv), and each barcode's status.
#
#   Rscript oracle.R identify_empty_drops.R INPUT_DIR OUT_DIR
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
inputs <- sort(list.files(args[2], pattern = "\\.tsv\\.gz$", full.names = TRUE))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
rows <- character()
for (path in inputs) {
  case <- sub("\\.tsv\\.gz$", "", basename(path))
  table <- read.table(gzfile(path), sep = "\t", col.names = c("barcode", "nf", "umi"), colClasses = c("character", "numeric", "numeric"))
  warned <- FALSE
  called <- withCallingHandlers(
    identify_empty_drops(table[, c("nf", "umi")]),
    warning = function(w) { if (grepl("more than one peak", conditionMessage(w))) warned <<- TRUE; invokeRestart("muffleWarning") }
  )
  # The bandwidth and the cut as identify_empty_drops computes them (it returns neither).
  fit <- ks::kdde(x = table$nf, deriv.order = 1)
  runs <- rle(fit$estimate > 0)
  cut <- fit$eval.points[sum(runs$lengths[1:2])]
  writeLines(called$cell_status, gzfile(file.path(args[3], paste0(case, ".status.gz"))))
  rows <- c(rows, sprintf("%s\t%.17g\t%.17g\t%d\t%s", case, fit$h, cut, length(runs$lengths), warned))
}
writeLines(c("case\th\tcut\truns\tone_peak_warning", rows), file.path(args[3], "summary.tsv"))
writeLines(c(sprintf("R %s", getRversion()), sprintf("ks %s", packageVersion("ks"))), file.path(args[3], "versions.txt"))
