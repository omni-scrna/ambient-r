#!/usr/bin/env Rscript
# DecontX ambient removal module.
#
# Reads a raw h5ad (one pool per dataset) + called cells TSV, subsets the raw
# data to called cells, runs decontX, and writes the corrected matrix as h5ad.

suppressPackageStartupMessages({
  library(decontX)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(data.table)
})

source("src/common/cli.R")
p <- arg_parser("DecontX ambient removal module")
p <- add_base_args(p)
p <- add_stage_args(p, "cleaning")
p <- add_argument(p, "--use_background", default = "true", help = "Whether to use raw matrix as background")
args <- parse_args(p, argv = strip_output_args())
args$use_background <- tolower(args$use_background) %in% c("true", "1")

# logging
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
  # reconstructed assay is named "X" instead of "counts" (decontX expects an
  # assay literally named "counts").
  sce_raw <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")
  sce_filt <- sce_raw[, cells_dt$cell_id]

  if (args$use_background) {
    message("  running decontX with raw matrix as background")
    dcx_res_sce <- decontX(sce_filt, background = sce_raw)
  } else {
    message("  running decontX in default mode")
    dcx_res_sce <- decontX(sce_filt)
  }
  corrected_mat <- round(assay(dcx_res_sce, "decontXcounts"))
  sce_out <- SingleCellExperiment(
    assays = list(counts = corrected_mat),
    rowData = rowData(sce_filt),
    colData = colData(sce_filt)
  )

  out_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
  # x_mapping: without it, anndataR leaves X empty and puts the assay in
  # layers["counts"] instead -- invisible to any plain scanpy.read_h5ad() consumer.
  write_h5ad(sce_out, out_path, x_mapping = "counts")
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
