#!/usr/bin/env Rscript
# Checks for prep.R on small synthetic pools. Run from the repo root:
#   Rscript tests/test_prep.R
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(data.table)
})
source("prep.R")  # defines the functions; main() only runs as a script

sce_of <- function(counts, rd, cd) {
  SingleCellExperiment(assays = list(counts = as(counts, "CsparseMatrix")), rowData = rd, colData = cd)
}

# --- orthologs: TP53/Trp53 and GAPDH/Gapdh are 1:1; one hashtag feature is dropped
genes <- c("GRCh38_ENSG00000141510", "GRCh38_ENSG00000111640",
           "GRCm39_ENSMUSG00000059552", "GRCm39_ENSMUSG00000057666", "HTO_1")
m <- matrix(1:20, nrow = 5, dimnames = list(genes, paste0("c", 1:4)))
mixed <- sce_of(m,
  DataFrame(species = c("human", "human", "mouse", "mouse", NA),
            feature_type = c(rep("Gene Expression", 4), "Antibody Capture"),
            Symbol = c("TP53", "GAPDH", "Trp53", "Gapdh", "HTO_1"), row.names = genes),
  DataFrame(pool_id = rep("p1", 4), row.names = colnames(m)))
res <- prepare(mixed, c("c2", "c3"), "orthologs")
stopifnot(
  setequal(rownames(res$sce), c("ENSG00000141510_ENSMUSG00000059552", "ENSG00000111640_ENSMUSG00000057666")),
  all(as.matrix(counts(res$sce))["ENSG00000141510_ENSMUSG00000059552", ] == m[1, ] + m[3, ]),
  all(as.matrix(counts(res$sce))["ENSG00000111640_ENSMUSG00000057666", ] == m[2, ] + m[4, ]),
  all(rowData(res$sce)$species == "ortholog"),
  identical(colnames(res$sce), colnames(m)),
  identical(res$cell_ids, c("c2", "c3"))
)

# --- gems: G1 was read under two probes, G2 under one
cells <- c("G1_BC001", "G1_BC002", "G2_BC001")
f <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, dimnames = list(c("g1", "g2"), cells))
flex <- sce_of(f,
  DataFrame(feature_type = rep("Gene Expression", 2), row.names = rownames(f)),
  DataFrame(pool_id = rep("p1", 3), probe_id = c("BC001", "BC002", "BC001"), row.names = cells))
res <- prepare(flex, c("G1_BC001", "G1_BC002"), "gems")
stopifnot(
  identical(colnames(res$sce), c("G1", "G2")),
  all(as.matrix(counts(res$sce))[, "G1"] == f[, 1] + f[, 2]),
  all(as.matrix(counts(res$sce))[, "G2"] == f[, 3]),
  !"probe_id" %in% colnames(colData(res$sce)),
  all(colData(res$sce)$pool_id == "p1"),
  identical(res$cell_ids, "G1")
)

# --- none: only the feature filter applies; a merge on the wrong pool is an error
res <- prepare(flex, cells, "none")
stopifnot(identical(dim(res$sce), dim(flex)), identical(res$cell_ids, cells))
stopifnot(inherits(try(prepare(flex, cells, "orthologs"), silent = TRUE), "try-error"))
stopifnot(inherits(try(prepare(mixed, "c1", "gems"), silent = TRUE), "try-error"))

# --- the script, end to end, then a method on its outputs
tmp <- tempfile("prep"); dir.create(tmp)
write_h5ad(flex, file.path(tmp, "raw.h5ad"), x_mapping = "counts")
fwrite(data.table(cell_id = c("G1_BC001", "G1_BC002")), file.path(tmp, "called.tsv"), sep = "\t")
run <- function(script, ...) {
  status <- system2("Rscript", c(script, "--output_dir", tmp, ...), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop(script, " failed")
}
for (merge in c("none", "gems")) {
  run("prep.R", "--name", merge, "--merge", merge, "--output", "prepared_h5ad=ignored",
      "--rawdata_h5ad", file.path(tmp, "raw.h5ad"), "--called_cells_tsv", file.path(tmp, "called.tsv"))
  run("passthrough.R", "--name", paste0("out_", merge),
      "--prepared_h5ad", file.path(tmp, paste0(merge, "_prepared.h5ad")),
      "--prepared_cells_tsv", file.path(tmp, paste0(merge, "_prepared_cells.tsv")))
}
none <- read_h5ad(file.path(tmp, "out_none.h5ad"), as = "SingleCellExperiment", x_mapping = "counts")
gems <- read_h5ad(file.path(tmp, "out_gems.h5ad"), as = "SingleCellExperiment", x_mapping = "counts")
stopifnot(
  identical(colnames(none), c("G1_BC001", "G1_BC002")),
  identical(colnames(gems), "G1"),
  all(rowSums(as.matrix(counts(none))) == as.matrix(counts(gems))[, "G1"])
)
cat("prep.R: all checks passed\n")
