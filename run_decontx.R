#!/usr/bin/env Rscript
# DecontX ambient removal module.
#
# Reads a raw h5ad (one pool per dataset) + called cells TSV, subsets the raw
# data to called cells, concatenated over orthologs for mixed species datasets,
# runs decontX, and writes the corrected matrix as h5ad.


suppressPackageStartupMessages({
  library(decontX)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
  library(data.table)
  library(babelgene)
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


merge_orthologs <- function(sce) {
  if (!"species" %in% colnames(rowData(sce))) {
    return(sce) 
  }
  message("  mixed-species dataset: merging orthologs before ambient removal ..")

  is_human <- rowData(sce)$species == "human"
  is_mouse <- rowData(sce)$species == "mouse"

  # strip the genome prefix to get ensembl_ids
  bare_geneids  <- sub("^[^_]+_", "", rownames(sce))
  human_geneids <- bare_geneids[is_human]
  mouse_geneids <- bare_geneids[is_mouse]

  # find othologs
  ortho_dt <- babelgene::orthologs(genes = human_geneids, species = "mouse", human = TRUE) %>%
    as.data.table() %>%
    setnames(c("ensembl", "symbol"), c("mouse_ensembl", "mouse_symbol")) %>%
    .[mouse_ensembl %in% mouse_geneids, ] %>%
    .[, n_human := .N, by = human_ensembl] %>%
    .[, n_mouse := .N, by = mouse_ensembl] %>%
    # drop any gene that has more than one match
    .[n_human == 1 & n_mouse == 1, .(human_ensembl, human_symbol, mouse_ensembl, mouse_symbol)]

  # merge orthologs
  human_rows <- setNames(which(is_human), human_geneids)
  mouse_rows <- setNames(which(is_mouse), mouse_geneids)

  human_rows_ortho <- human_rows[ortho_dt$human_ensembl]
  mouse_rows_ortho <- mouse_rows[ortho_dt$mouse_ensembl]

  merged_counts <- counts(sce)[human_rows_ortho, , drop = FALSE] + counts(sce)[mouse_rows_ortho, , drop = FALSE]
  rownames(merged_counts) <- paste(ortho_dt$human_ensembl, ortho_dt$mouse_ensembl, sep = '_')

  # sort out rowData
  merged_rd <- rowData(sce)[human_rows_ortho, , drop = FALSE]
  rownames(merged_rd) <- rownames(merged_counts)
  merged_rd$species <- "ortholog"
  if ("ID" %in% colnames(merged_rd)) merged_rd$ID <- rownames(merged_rd)
  if ("Symbol" %in% colnames(merged_rd)) merged_rd$Symbol <- paste(ortho_dt$human_symbol, ortho_dt$mouse_symbol, sep = '_')

  # return sce with merged counts (the output only ever contains these merged rows, every
  # non-orthologous human/mouse gene is dropped)
  SingleCellExperiment(
    assays = list(counts = merged_counts),
    rowData = merged_rd,
    colData = colData(sce)
  )
}

merge_gems <- function(sce) {
  if (!"probe_id" %in% colnames(colData(sce))) {
    return(sce)
  }

  # column name is "{gem_barcode}_{probe_id}"
  gem_id <- sub("_[^_]+$", "", colnames(sce))
  if (!any(duplicated(gem_id))) {
    message("  flex dataset: no GEM barcode appears under more than one probe_id -- nothing to merge")
    # still drop to GEM-level colnames for consistency with the
    colnames(sce) <- gem_id
    colData(sce)$probe_id <- NULL
    return(sce)
  }
  message("  flex dataset: merging raw counts for columns sharing the same GEM barcode before ambient removal ..")

  # sum counts for every column sharing a GEM barcode into one column
  ind <- Matrix::fac2sparse(gem_id)
  merged_counts <- counts(sce) %*% Matrix::t(ind)
  colnames(merged_counts) <- rownames(ind)

  # one colData row per distinct GEM
  first_idx <- match(colnames(merged_counts), gem_id)
  merged_cd <- colData(sce)[first_idx, , drop = FALSE]
  rownames(merged_cd) <- colnames(merged_counts)
  merged_cd$probe_id <- NULL

  message(sprintf("  merged %d raw GEM x probe_id columns into %d GEM barcodes",
    ncol(sce), ncol(merged_counts)))

  SingleCellExperiment(
    assays = list(counts = merged_counts),
    rowData = rowData(sce),
    colData = merged_cd
  )
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading input files ..")
  cells_dt <- fread(args$called_cells_tsv)
  sce_raw <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message("  keeping only RNA (Gene Expression) features ..")
  sce_raw <- sce_raw[rowData(sce_raw)$feature_type == "Gene Expression", ]
 
  # merge othologs
  sce_raw <- merge_orthologs(sce_raw)

  # merge counts for each GEM in flex data -- check before merge_gems() runs,
  # since it strips colData$probe_id either way
  is_flex <- "probe_id" %in% colnames(colData(sce_raw))
  sce_raw <- merge_gems(sce_raw)

  if (is_flex) {
    # cell_calling writes cell_id as "{gem_barcode}_{probe_id}"
    # strip to GEM barcode to match sce_raw's new GEM-level colnames.
    called_ids <- unique(sub("_[^_]+$", "", cells_dt$cell_id))
  } else {
    called_ids <- cells_dt$cell_id
  }
  sce_filt <- sce_raw[, called_ids]

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
