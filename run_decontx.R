#!/usr/bin/env Rscript
# DecontX ambient removal module.
#
# Reads raw 10x H5 + filtered h5ad matrices, runs decontX with the raw matrix
# as background or in the default mode, and writes a corrected h5ad matrix.

suppressPackageStartupMessages({
  library(decontX)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(Matrix)
  library(anndataR)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))


main <- function() {
  args <- parse_args_checked(description = "DecontX ambient removal module")
  message(sprintf("Full command: %s", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in names(args)) message(sprintf("  %s: %s", k, args[[k]]))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading filtered matrix")
  sce_filt <- read_h5ad(args$filtered_h5ad, as = "SingleCellExperiment")

  if (args$use_background) {
    message("  reading raw matrix (background)")
    sce_raw <- read10xCounts(args$rawdata_h5, type = "HDF5")

    message("  running decontX with background")
    sce_result <- decontX(sce_filt, background = sce_raw)
  } else {
    message("  running decontX without background (default)")
    sce_result <- decontX(sce_filt)
  }

  corrected <- assay(sce_result, "decontXcounts")
  corrected <- round(corrected)

  sce_out <- SingleCellExperiment(
    assays = list(X = corrected),
    rowData = rowData(sce_filt),
    colData = colData(sce_filt)
  )

  out_path <- file.path(args$output_dir, paste0(args$name, "_corrected.h5ad"))
  write_h5ad(sce_out, out_path)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
