#!/usr/bin/env Rscript
# DecontX ambient removal module.
#
# Reads the PREP outputs of one pool (count matrix + called cells; raw, or
# merged over orthologs / GEMs by PREP), subsets the matrix to the called
# cells, runs decontX, and writes the corrected matrix as h5ad.


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
p <- add_stage_args(p, "AMBIENT")
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
  cells_dt <- fread(args$prepared_cells_tsv)
  sce_raw <- read_h5ad(args$prepared_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message("  keeping only RNA (Gene Expression) features ..")
  sce_raw <- sce_raw[rowData(sce_raw)$feature_type == "Gene Expression", ]
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
  write_h5ad(sce_out, out_path, x_mapping = "counts")
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
