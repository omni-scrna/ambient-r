#!/usr/bin/env Rscript
# SoupX ambient removal module.
#
# Reads raw 10x H5 + filtered h5ad matrices, runs SoupX to estimate and remove
# ambient RNA contamination, and writes a corrected h5ad matrix.

suppressPackageStartupMessages({
  library(SoupX)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(scran)
  library(scuttle)
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
  args <- parse_args_checked(description = "SoupX ambient removal module")
  message(sprintf("Full command: %s", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in names(args)) message(sprintf("  %s: %s", k, args[[k]]))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading raw matrix")
  sce_raw <- read10xCounts(args$rawdata_h5, type = "HDF5")

  message("  reading filtered matrix")
  sce_filt <- read_h5ad(args$filtered_h5ad) |> to_SingleCellExperiment()
  message(sprintf("  filtered: %d genes x %d cells", nrow(sce_filt), ncol(sce_filt)))

  # quick clustering for SoupX (requires cluster labels)
  message("  clustering for SoupX")
  sce_pp <- logNormCounts(sce_filt)
  dec <- modelGeneVar(sce_pp)
  hvgs <- getTopHVGs(dec, n = 2000)
  sce_pp <- runPCA(sce_pp, subset_row = hvgs)
  g <- buildSNNGraph(sce_pp, use.dimred = "PCA")
  clusters <- igraph::cluster_louvain(g)$membership
  names(clusters) <- colnames(sce_filt)

  # build SoupX channel
  message("  running SoupX")
  raw_mat <- counts(sce_raw)
  filt_mat <- assay(sce_filt, "X")
  rownames(raw_mat) <- rowData(sce_raw)$Symbol
  rownames(filt_mat) <- rowData(sce_filt)$Symbol
  colnames(raw_mat) <- colData(sce_raw)$Barcode
  colnames(filt_mat) <- colnames(sce_filt)

  sc <- SoupChannel(raw_mat, filt_mat, calcSoupProfile = FALSE)
  soup_prof <- data.frame(
    row.names = rownames(filt_mat),
    est = rowSums(raw_mat) / sum(raw_mat),
    counts = rowSums(raw_mat)
  )
  sc <- setSoupProfile(sc, soup_prof)
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  message(sprintf("  estimated contamination (rho): %.4f", sc$metaData$rho[1]))

  corrected <- adjustCounts(sc, roundToInt = TRUE)

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
