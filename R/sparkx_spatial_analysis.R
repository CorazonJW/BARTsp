#' Run SPARKX for Spatial Gene Expression Analysis
#'
#' This function runs SPARKX to identify spatially variable genes from spatial transcriptomics data.
#'
#' @param object A list containing \code{expression_matrix} (gene expression data) 
#' and \code{spatial_coordinates} (spatial locations).
#' @param numCores Number of CPU cores to use for computation.
#'
#' @return A SPARKX result object containing spatially variable gene analysis results.
#'
#' @importFrom SPARK sparkx
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @export
run_SPARKX <- function(object, numCores) {
  # Input validation
  if (!is.list(object)) {
    stop("Error: object must be a list.")
  }
  
  if (!all(c("expression_matrix", "spatial_coordinates") %in% names(object))) {
    stop("Error: object must contain 'expression_matrix' and 'spatial_coordinates'.")
  }
  
  if (!is.matrix(object$expression_matrix) && !is.data.frame(object$expression_matrix)) {
    stop("Error: expression_matrix must be a matrix or data frame.")
  }
  
  if (!is.data.frame(object$spatial_coordinates)) {
    stop("Error: spatial_coordinates must be a data frame.")
  }
  
  if (!all(c("barcodes", "x", "y") %in% colnames(object$spatial_coordinates))) {
    stop("Error: spatial_coordinates must contain 'barcodes', 'x', and 'y' columns.")
  }
  
  if (!is.numeric(numCores) || numCores < 1) {
    stop("Error: numCores must be a positive integer.")
  }

  # Input
  sp_count <- object$expression_matrix
  coords_df <- object$spatial_coordinates

  # Ensure spatial coordinates match expression matrix
  matching_cells <- intersect(colnames(sp_count), coords_df$barcodes)
  if (length(matching_cells) == 0) {
    stop("Error: No matching cells found between 'expression_matrix' and 'spatial_coordinates'.")
  }
  
  # Subset to matching cells
  sp_count <- sp_count[, matching_cells, drop = FALSE]
  coords_df <- coords_df[coords_df$barcodes %in% matching_cells, , drop = FALSE]

  # Ensure coordinates match column order
  coords_df <- coords_df[match(colnames(sp_count), coords_df$barcodes), ]

  # Check for missing or invalid coordinates
  if (any(is.na(coords_df$x)) || any(is.na(coords_df$y))) {
    stop("Error: Spatial coordinates contain missing values.")
  }
  
  if (any(!is.numeric(coords_df$x)) || any(!is.numeric(coords_df$y))) {
    stop("Error: Spatial coordinates must be numeric.")
  }

  # Rename expression matrix
  colnames(sp_count) <- paste0(coords_df$x, "x", coords_df$y)

  # Prepare location matrix
  location <- tryCatch({
    info <- cbind.data.frame(
      x = as.numeric(sapply(strsplit(colnames(sp_count), split = "x"), "[", 1)),
      y = as.numeric(sapply(strsplit(colnames(sp_count), split = "x"), "[", 2))
    )
    rownames(info) <- colnames(sp_count)
    as.matrix(info)
  }, error = function(e) {
    stop("Error processing spatial coordinates: ", e$message)
  })

  # Remove mitochondrial genes safely
  mt_idx <- grep("^mt-", rownames(sp_count), ignore.case = TRUE)
  if (length(mt_idx) > 0) {
    sp_count <- sp_count[-mt_idx, , drop = FALSE]
  }

  # Ensure expression matrix is not empty
  if (nrow(sp_count) == 0) {
    stop("Error: No genes left after filtering mitochondrial genes.")
  }
  if (ncol(sp_count) == 0) {
    stop("Error: No cells in the expression matrix after processing.")
  }

  # Check for zero variance genes
  zero_var_genes <- apply(sp_count, 1, function(x) var(x) == 0)
  if (any(zero_var_genes)) {
    warning("Removing ", sum(zero_var_genes), " genes with zero variance.")
    sp_count <- sp_count[!zero_var_genes, , drop = FALSE]
  }

  # Run SPARKX with error handling
  sparkX_result <- tryCatch({
    SPARK::sparkx(sp_count, location, numCores = numCores, option = "mixture")
  }, error = function(e) {
    stop("Error during SPARKX execution: ", e$message)
  })

  return(sparkX_result)
}

#' Extract Differentially Expressed Genes from SPARKX Results
#'
#' This function extracts significant spatially variable genes from a SPARKX analysis.
#'
#' @param spark_obj A SPARKX result object.
#' @param cutoff Adjusted p-value threshold for significance.
#'
#' @return A character vector of significant spatially variable genes.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @export
get_sparkx_DEGs <- function(spark_obj, cutoff) {
  # Input validation
  if (!is.list(spark_obj)) {
    stop("Error: spark_obj must be a list.")
  }
  
  if (!"res_mtest" %in% names(spark_obj)) {
    stop("Error: spark_obj must contain 'res_mtest' element.")
  }
  
  if (!is.numeric(cutoff) || cutoff < 0 || cutoff > 1) {
    stop("Error: cutoff must be a numeric value between 0 and 1.")
  }
  
  if (!"adjustedPval" %in% colnames(spark_obj$res_mtest)) {
    stop("Error: spark_obj$res_mtest must contain 'adjustedPval' column.")
  }

  # Extract significant genes
  sig_results <- tryCatch({
    spark_obj$res_mtest %>% dplyr::filter(adjustedPval < cutoff)
  }, error = function(e) {
    stop("Error filtering results: ", e$message)
  })
  
  significant_genes <- rownames(sig_results)
  
  if (length(significant_genes) == 0) {
    warning("No significant genes found at the specified cutoff.")
  }

  return(list(significant_features = significant_genes))
} 