suppressPackageStartupMessages(library(optparse))

build_parser <- function(description = "Ambient removal module") {
  option_list <- list(
    make_option("--output_dir", type = "character",
                help = "Output directory for results"),
    make_option("--name", type = "character",
                help = "Module name/identifier"),
    make_option("--rawdata.h5", type = "character",
                help = "Raw 10x H5 count matrix (unfiltered)"),
    make_option("--filtered.h5ad", type = "character",
                help = "Filtered h5ad count matrix (cells only)"),
    make_option("--use_background", type = "character", default = "true",
                help = "Use raw matrix as background [default: true]")
  )
  OptionParser(option_list = option_list, description = description)
}

parse_args_checked <- function(description = "Ambient removal module",
                               require_raw = TRUE) {
  parser <- build_parser(description)
  raw <- parse_args(parser)

  args <- list(
    output_dir    = raw$output_dir,
    name          = raw$name,
    rawdata_h5    = raw[["rawdata.h5"]],
    filtered_h5ad = raw[["filtered.h5ad"]],
    use_background = tolower(raw$use_background) %in% c("true", "yes", "1")
  )

  required <- c("output_dir", "name", "filtered_h5ad")
  if (require_raw) required <- c(required, "rawdata_h5")
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  args
}
