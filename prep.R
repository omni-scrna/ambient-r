#!/usr/bin/env Rscript
# PREP: the count matrix every ambient removal method of one lineage runs on.
#
# Reads the raw h5ad of one pool + the called cells TSV, and writes both again
# in the id space the methods will see. `--merge` picks the variant:
#
#   none       the raw matrix and the called cells, unchanged
#   orthologs  mixed-species pools: human and mouse 1:1 orthologs are summed
#              into one feature; every other gene is dropped
#   gems       Flex pools: the GEM x probe_id columns of one GEM barcode are
#              summed into one column, and cell ids lose their "_{probe_id}"
#
# The plan gates the merges with `requires` on the dataset's `multiplexing`
# label, so a merge that does not apply to its input is an error here.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(data.table)
})

MERGES <- c("none", "orthologs", "gems")

merge_orthologs <- function(sce) {
  if (!"species" %in% colnames(rowData(sce))) {
    stop("--merge orthologs needs rowData 'species' (a mixed-species pool)")
  }
  is_human <- rowData(sce)$species == "human"
  is_mouse <- rowData(sce)$species == "mouse"

  # strip the genome prefix to get ensembl ids
  bare_geneids <- sub("^[^_]+_", "", rownames(sce))
  human_rows <- setNames(which(is_human), bare_geneids[is_human])
  mouse_rows <- setNames(which(is_mouse), bare_geneids[is_mouse])

  ortho_dt <- as.data.table(
    babelgene::orthologs(genes = names(human_rows), species = "mouse", human = TRUE)
  )
  setnames(ortho_dt, c("ensembl", "symbol"), c("mouse_ensembl", "mouse_symbol"))
  ortho_dt <- ortho_dt[mouse_ensembl %in% names(mouse_rows)]
  ortho_dt[, n_human := .N, by = human_ensembl]
  ortho_dt[, n_mouse := .N, by = mouse_ensembl]
  # 1:1 pairs only
  ortho_dt <- ortho_dt[n_human == 1 & n_mouse == 1]
  if (nrow(ortho_dt) == 0) stop("no 1:1 human-mouse orthologs found in this matrix")

  human_rows_ortho <- human_rows[ortho_dt$human_ensembl]
  mouse_rows_ortho <- mouse_rows[ortho_dt$mouse_ensembl]
  merged_counts <- counts(sce)[human_rows_ortho, , drop = FALSE] +
    counts(sce)[mouse_rows_ortho, , drop = FALSE]
  rownames(merged_counts) <- paste(ortho_dt$human_ensembl, ortho_dt$mouse_ensembl, sep = "_")

  merged_rd <- rowData(sce)[human_rows_ortho, , drop = FALSE]
  rownames(merged_rd) <- rownames(merged_counts)
  merged_rd$species <- "ortholog"
  if ("ID" %in% colnames(merged_rd)) merged_rd$ID <- rownames(merged_rd)
  if ("Symbol" %in% colnames(merged_rd)) {
    merged_rd$Symbol <- paste(ortho_dt$human_symbol, ortho_dt$mouse_symbol, sep = "_")
  }
  message(sprintf("  merged %d features into %d ortholog pairs", nrow(sce), nrow(merged_counts)))

  SingleCellExperiment(
    assays = list(counts = merged_counts),
    rowData = merged_rd,
    colData = colData(sce)
  )
}

# "{gem_barcode}_{probe_id}" -> "{gem_barcode}"
gem_of <- function(cell_id) sub("_[^_]+$", "", cell_id)

merge_gems <- function(sce) {
  if (!"probe_id" %in% colnames(colData(sce))) {
    stop("--merge gems needs colData 'probe_id' (a Flex pool)")
  }
  gem_id <- gem_of(colnames(sce))
  # sum the columns sharing a GEM barcode into one
  ind <- Matrix::fac2sparse(gem_id)
  merged_counts <- counts(sce) %*% Matrix::t(ind)
  colnames(merged_counts) <- rownames(ind)

  # one colData row per GEM; probe_id no longer describes a column
  merged_cd <- colData(sce)[match(colnames(merged_counts), gem_id), , drop = FALSE]
  rownames(merged_cd) <- colnames(merged_counts)
  merged_cd$probe_id <- NULL
  message(sprintf("  merged %d GEM x probe_id columns into %d GEM barcodes",
                  ncol(sce), ncol(merged_counts)))

  SingleCellExperiment(
    assays = list(counts = merged_counts),
    rowData = rowData(sce),
    colData = merged_cd
  )
}

# Returns list(sce, cell_ids): the matrix and the called cells, in the same id space.
prepare <- function(sce, cell_ids, merge) {
  sce <- sce[rowData(sce)$feature_type == "Gene Expression", ]
  if (merge == "orthologs") {
    sce <- merge_orthologs(sce)
  } else if (merge == "gems") {
    sce <- merge_gems(sce)
    cell_ids <- unique(gem_of(cell_ids))
  }
  missing <- setdiff(cell_ids, colnames(sce))
  if (length(missing) > 0) {
    stop(sprintf("%d called cell(s) are not in the prepared matrix, e.g. %s",
                 length(missing), missing[1]))
  }
  list(sce = sce, cell_ids = cell_ids)
}

main <- function() {
  source("src/common/cli.R")
  p <- arg_parser("PREP: raw or merged counts for the ambient removal methods")
  p <- add_base_args(p)
  p <- add_stage_args(p, "PREP")
  p <- add_argument(p, "--merge", default = "none",
                    help = "none | orthologs (mixed-species) | gems (Flex)")
  args <- parse_args(p, argv = strip_output_args())
  if (!args$merge %in% MERGES) {
    stop(sprintf("--merge must be one of: %s", paste(MERGES, collapse = ", ")))
  }
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  out_h5ad <- file.path(args$output_dir, paste0(args$name, "_prepared.h5ad"))
  out_tsv <- file.path(args$output_dir, paste0(args$name, "_prepared_cells.tsv"))

  if (args$merge == "none") {
    # Byte-identical to the input: a hard link where the filesystem allows it.
    unlink(c(out_h5ad, out_tsv))
    if (!suppressWarnings(file.link(args$rawdata_h5ad, out_h5ad))) {
      file.copy(args$rawdata_h5ad, out_h5ad)
    }
    file.copy(args$called_cells_tsv, out_tsv)
  } else {
    message("  reading input files ..")
    sce_raw <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")
    res <- prepare(sce_raw, fread(args$called_cells_tsv)$cell_id, args$merge)
    write_h5ad(res$sce, out_h5ad, x_mapping = "counts")
    fwrite(data.table(cell_id = res$cell_ids), out_tsv, sep = "\t", quote = FALSE)
  }
  message(sprintf("  wrote: %s, %s", out_h5ad, out_tsv))
}

if (sys.nframe() == 0L) {
  main()
}
