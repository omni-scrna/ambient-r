#!/usr/bin/env Rscript
# scCDC ambient removal module.
#
# Reads a filtered h5ad matrix, runs scCDC contamination detection and
# correction via Seurat, and writes a corrected h5ad matrix.

suppressPackageStartupMessages({
  library(Seurat)
  library(scCDC)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(Matrix)
  library(anndataR)
})

source("src/common/cli.R")
p <- arg_parser("scCDC ambient removal module")
p <- add_base_args(p)
p <- add_stage_args(p, "cleaning")
p <- add_argument(p, "--filtered_h5ad", help = "Filtered h5ad matrix")
args <- parse_args(p)

message(sprintf("Full command: %s", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
for (k in names(args)) message(sprintf("  %s: %s", k, args[[k]]))


main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading filtered matrix ..")
  sce_filt <- read_h5ad(args$filtered_h5ad) |> to_SingleCellExperiment()
  mat <- assay(sce_filt, "X")
  rownames(mat) <- rowData(sce_filt)$Symbol
  colnames(mat) <- colnames(sce_filt)
  message(sprintf("  filtered: %d genes x %d cells", nrow(mat), ncol(mat)))

  # filter genes with zero cells
  keep_genes <- rowSums(mat > 0) >= 1
  mat <- mat[keep_genes, ]
  message(sprintf("  after gene filter (min_cells=1): %d genes", nrow(mat)))

  # Seurat preprocessing required by scCDC
  message("  creating Seurat object and preprocessing ..")
  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                              scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                     nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj,
                       features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

  message("  running scCDC contamination detection ..")
  GCGs <- ContaminationDetection(seurat_obj)
  message(sprintf("  detected %d contamination genes", nrow(GCGs)))

  message("  running scCDC contamination correction ..")
  seurat_corrected <- ContaminationCorrection(seurat_obj, rownames(GCGs))

  DefaultAssay(seurat_corrected) <- "Corrected"
  corrected <- GetAssayData(seurat_corrected, layer = "counts")
  message(sprintf("  corrected matrix: %d genes x %d cells",
              nrow(corrected), ncol(corrected)))

  # map back to original gene metadata
  gene_idx <- match(rownames(corrected), rowData(sce_filt)$Symbol)
  gene_ids <- rowData(sce_filt)$ID[gene_idx]
  gene_symbols <- rowData(sce_filt)$Symbol[gene_idx]

  sce_out <- SingleCellExperiment(
    assays = list(X = corrected),
    rowData = DataFrame(ID = gene_ids, Symbol = gene_symbols),
    colData = colData(sce_filt)[match(colnames(corrected), colnames(sce_filt)), ]
  )

  out_path <- file.path(args$output_dir, paste0(args$name, "_corrected.h5ad"))
  adata <- from_SingleCellExperiment(sce_out)
  write_h5ad(adata, out_path)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
