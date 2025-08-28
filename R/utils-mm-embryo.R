# Internal checksum + helpers for large dataset

# MD5 for mm_embryo_coprofiling.rda (update if you replace the asset)
.mm_embryo_expected_md5 <- "37a38699e01f27b0b39e296cfebe6b29"

.verify_md5 <- function(path, expected) {
  if (is.null(expected)) return(invisible(TRUE))
  has <- as.character(tools::md5sum(path))
  if (!identical(has, expected)) {
    stop(
      "Checksum mismatch for '", basename(path), "'.\n",
      "Expected: ", expected, "\nGot:      ", has, "\n",
      "Delete the cached file and try again:\n",
      "BARTsp::bartsp_clear_cache(); BARTsp::load_mm_embryo()"
    )
  }
  invisible(TRUE)
}

#' Directory where large BARTsp datasets are cached
#' @export
bartsp_data_dir <- function() rappdirs::user_data_dir("BARTsp")

#' Clear cached large datasets
#'
#' Deletes files under BARTsp's user cache directory.
#' @param ask Prompt before deleting?
#' @export
bartsp_clear_cache <- function(ask = interactive()) {
  d <- bartsp_data_dir()
  if (!dir.exists(d)) return(invisible(TRUE))
  if (ask && !utils::askYesNo(paste0("Delete cache at ", d, "?"))) return(invisible(FALSE))
  unlink(d, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

#' Check if full mm_embryo_coprofiling is already cached
#' @export
has_mm_embryo <- function() {
  file.exists(file.path(bartsp_data_dir(), "mm_embryo_coprofiling.rda"))
}
