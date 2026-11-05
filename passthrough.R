#!/usr/bin/env Rscript
# Passthrough "no ambient removal" module.
#
# Baseline for the cleaning stage: subsets the raw h5ad to the called-cell
# barcodes with no correction applied, for comparison against decontx/cellbender.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(data.table)
})

source("src/common/cli.R")
p <- arg_parser("Passthrough (no ambient removal) module")
p <- add_base_args(p)
p <- add_stage_args(p, "cleaning")
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading input files ..")
  cells_dt <- fread(args$called_cells_tsv)
  # x_mapping: must match the x_mapping="counts" used on write, or the
  # reconstructed assay is named "X" instead of "counts".
  sce_raw <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message(sprintf("  subsetting to %d called cells (no correction applied)", nrow(cells_dt)))
  sce_filt <- sce_raw[, cells_dt$cell_id]

  out_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
  # x_mapping: without it, anndataR leaves X empty and puts the assay in
  # layers["counts"] instead -- invisible to any plain scanpy.read_h5ad() consumer.
  write_h5ad(sce_filt, out_path, x_mapping = "counts")
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
