#' Load the full mm_embryo_coprofiling dataset (download-once cache)
#'
#' This helper downloads a large `.rda` from a GitHub Release (or a direct URL),
#' caches it in the user's data directory, and loads the object into the current session.
#'
#' @param dest Cache directory; default is a per-user cache.
#' @param repo GitHub "owner/repo" string used by piggyback.
#' @param tag  GitHub Release tag that contains the asset.
#' @param asset Filename of the asset on the Release.
#' @param download_url Optional direct URL to the .rda (if you don't want piggyback).
#' @param envir Environment to load the object into. Default: `.GlobalEnv`.
#' @return (Invisibly) the loaded `mm_embryo_coprofiling` object.
#' @export
#' @examples
#' \dontrun{
#'   mm <- load_mm_embryo()
#'   str(mm)
#' }
load_mm_embryo <- function(
  dest = rappdirs::user_data_dir("BARTsp"),
  repo = "CorazonJW/BARTsp",
  tag  = "v0.1-data",
  asset = "mm_embryo_coprofiling.rda",
  download_url = NULL,
  envir = .GlobalEnv
) {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  local_file <- file.path(dest, asset)

  if (!file.exists(local_file)) {
    message("mm_embryo_coprofiling not cached; downloading...")

    if (!is.null(download_url)) {
      utils::download.file(download_url, local_file, mode = "wb", quiet = TRUE)
    } else if (requireNamespace("piggyback", quietly = TRUE)) {
      # Use piggyback to fetch from GitHub Releases
      piggyback::pb_download(
        file = asset, repo = repo, tag = tag, dest = dest, overwrite = TRUE
      )
    } else {
      stop(
        "Need to download the dataset, but {piggyback} isn't installed and no direct URL was provided.\n",
        "Install piggyback: install.packages('piggyback')\n",
        "Or call load_mm_embryo(download_url = 'https://.../mm_embryo_coprofiling.rda')"
      )
    }
  }

  # Load and verify object name
  obj_names <- load(local_file, envir = envir)
  if (!"mm_embryo_coprofiling" %in% obj_names) {
    stop(
      "Downloaded file did not contain an object named 'mm_embryo_coprofiling'. ",
      "It contained: ", paste(obj_names, collapse = ", ")
    )
  }

  invisible(get("mm_embryo_coprofiling", envir = envir))
}
