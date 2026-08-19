#!/usr/bin/env Rscript
# DecontX ambient removal module.
#
# Reads raw h5ad (with pool_id) + called cells TSV (with pool_id, sample_id),
# subsets raw data to called cells, runs decontX per pool, and writes DATA-compatible outputs.

suppressPackageStartupMessages({
  library(decontX)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(yaml)
  library(data.table)
})

source("src/common/cli.R")
p <- arg_parser("DecontX ambient removal module")
p <- add_base_args(p)
p <- add_stage_args(p, "cleaning")
p <- add_argument(p, "--use_background", default = "true", help = "Whether to use raw matrix as background")
args <- parse_args(p)
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
  sce_raw <- read_h5ad(args$rawdata_raw_h5ad, as = "SingleCellExperiment")

  pools = cells_dt$pool_id |> unique()

  res_ls = lapply(pools, function(pool) {
    message(sprintf("  --- processing pool: %s ---", pool))
    pool_cells <- cells_dt[pool_id == pool, cell_id]
    pool_filt <- sce_raw[, pool_cells]

    if (args$use_background) {
      pool_raw <- sce_raw[, colData(sce_raw)$pool_id == pool]
      message("  running decontX with raw matrix as background")
      dcx_res_sce <- decontX(pool_filt, background = pool_raw)
    } else {
      message("  running decontX in default mode")
      dcx_res_sce <- decontX(pool_filt)
    }
    corrected_mat <- round(assay(dcx_res_sce, "decontXcounts"))
    pool_out <- SingleCellExperiment(
      assays = list(X = corrected_mat),
      rowData = rowData(pool_filt),
      colData = colData(pool_filt)
    )
    
    return(pool_out)
  })


  message("  making sce with corrected counts for all pools ..")
  sce_combined <- do.call(cbind, res_ls)
  assayNames(sce_combined)[assayNames(sce_combined) == "X"] <- "counts"

  # write outputs
  h5ad_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
  write_h5ad(sce_combined, h5ad_path)
  message(sprintf("  wrote: %s", h5ad_path))

  truth_values <- as.character(colData(sce_combined)$sample_id)

  clusters_truth_path <- file.path(
    args$output_dir, paste0(args$name, ".clusters_truth.tsv")
  )
  fwrite(
    data.table(
      cell_id = colnames(sce_combined),
      truths = truth_values,
      truths_cl = NA_character_
    ),
    clusters_truth_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message(sprintf("  wrote: %s", clusters_truth_path))

  num_clusters_path <- file.path(
    args$output_dir, paste0(args$name, ".clusters_truth_num.txt")
  )
  writeLines(as.character(length(unique(truth_values))), con = num_clusters_path)
  message(sprintf("  wrote: %s", num_clusters_path))

  properties_path <- file.path(
    args$output_dir, paste0(args$name, "_properties.yaml")
  )
  write_yaml(
    list(
      batch_var = "pool_id",
      sample_var = "sample_id",
      labels_var = "sample_id"
    ),
    properties_path
  )
  message(sprintf("  wrote: %s", properties_path))
  
  
}

if (sys.nframe() == 0L) {
  main()
}
