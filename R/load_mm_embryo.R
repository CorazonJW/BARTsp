#' Load the full mm_embryo_coprofiling dataset (download-once cache)
#'
#' Downloads/caches a large `.rda` from a GitHub Release (or a direct URL),
#' then loads the object into the current session. If the file's internal
#' object name is not `mm_embryo_coprofiling` (e.g., it's `subset`), this
#' function will load it and assign it under the canonical name.
#'
#' @param dest Cache directory (default: per-user cache).
#' @param repo GitHub "owner/repo" used by \pkg{piggyback}.
#' @param tag  GitHub Release tag that contains the asset.
#' @param asset Filename of the asset on the Release.
#' @param download_url Optional direct URL to the .rda (skip \pkg{piggyback}).
#' @param envir Environment to load the object into (default: `.GlobalEnv`).
#' @param object_name Canonical name to assign in the session.
#' @return (Invisibly) the loaded object (also available as `object_name` in `envir`).
#' @export
load_mm_embryo <- function(
  dest = bartsp_data_dir(),
  repo = "CorazonJW/BARTsp",
  tag  = "v0.1-data",
  asset = "mm_embryo_coprofiling.rda",
  download_url = NULL,
  envir = .GlobalEnv,
  object_name = "mm_embryo_coprofiling"
) {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  local_file <- file.path(dest, asset)

  # Download if needed
  if (!file.exists(local_file)) {
    message(object_name, " not cached; downloading...")
    if (!is.null(download_url)) {
      utils::download.file(download_url, local_file, mode = "wb", quiet = TRUE)
    } else if (requireNamespace("piggyback", quietly = TRUE)) {
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

  # Load into a temporary environment to inspect object names
  tmp <- new.env(parent = emptyenv())
  obj_names <- load(local_file, envir = tmp)

  pick <- NULL
  if (object_name %in% obj_names) {
    pick <- object_name
  } else if ("subset" %in% obj_names) {
    pick <- "subset"
  } else if (length(obj_names) == 1L) {
    pick <- obj_names[[1L]]
  }

  if (is.null(pick)) {
    stop(
      "Downloaded file did not contain '", object_name, "'. ",
      "It contained: ", paste(obj_names, collapse = ", "), "."
    )
  }

  # Assign under the canonical name in target environment
  assign(object_name, get(pick, envir = tmp), envir = envir)
  invisible(get(object_name, envir = envir))
}
