#' Prepare input
#'
#' This function prepares input for following analysis.
#'
#' @param expression_matrix A sparse matrix of class \code{dgCMatrix} containing gene expression data.
#' @param cell_metadata A data frame containing metadata for each cell, including cell types.
#' @param feature_name A data frame containing feature metadata.
#' @param spatial_coordinates A data frame with 'x' and 'y' spatial coordinates for cells.
#' @param cell_types A character vector of cell types to select.
#' 
#' @return A list containing:
#' \item{expression_matrix}{Subsetted expression matrix for selected cells}
#' \item{cell_metadata}{Metadata for selected cells}
#' \item{feature_metadata}{Feature (gene) metadata}
#' \item{spatial_coordinates}{Spatial coordinates of selected cells}
#' 
#' @export
prepare_input <- function(expression_matrix, cell_metadata, feature_name, spatial_coordinates, cell_types) {
  # Input validation
  if (missing(expression_matrix) || missing(cell_metadata) || missing(feature_name) || 
      missing(spatial_coordinates) || missing(cell_types)) {
    stop("Error: All arguments are required.")
  }
  
  if (!inherits(expression_matrix, "dgCMatrix")) {
    stop("Error: expression_matrix must be a sparse matrix of class 'dgCMatrix'.")
  }
  
  if (!is.data.frame(cell_metadata)) {
    stop("Error: cell_metadata must be a data frame.")
  }
  
  if (!("cell_type" %in% colnames(cell_metadata))) {
    stop("Error: cell_metadata must contain a 'cell_type' column.")
  }
  
  if (!is.data.frame(feature_name)) {
    stop("Error: feature_metadata must be a data frame.")
  }

  if (!is.data.frame(spatial_coordinates)) {
    stop("Error: spatial_coordinates must be a data frame.")
  }
  
  if (!all(c("x", "y") %in% colnames(spatial_coordinates))) {
    stop("Error: spatial_coordinates must contain 'x' and 'y' columns.")
  }
  
  if (!is.character(cell_types) || length(cell_types) == 0) {
    stop("Error: cell_types must be a non-empty character vector.")
  }
  
  # Check for missing values
  if (any(is.na(cell_metadata$cell_type))) {
    stop("Error: cell_metadata contains missing values in the 'cell_type' column.")
  }
  
  if (any(is.na(spatial_coordinates$x)) || any(is.na(spatial_coordinates$y))) {
    stop("Error: spatial_coordinates contains missing values in 'x' or 'y' columns.")
  }
  
  # Check for non-numeric coordinates
  if (!is.numeric(spatial_coordinates$x) || !is.numeric(spatial_coordinates$y)) {
    stop("Error: spatial_coordinates 'x' and 'y' must be numeric.")
  }
  
  # Check for infinite values
  if (any(is.infinite(spatial_coordinates$x)) || any(is.infinite(spatial_coordinates$y))) {
    stop("Error: spatial_coordinates contains infinite values.")
  }

  # Extract cell barcodes matching the cell types
  cells <- tryCatch({
    rownames(cell_metadata)[cell_metadata$cell_type %in% cell_types]
  }, error = function(e) {
    stop("Error extracting cells: ", e$message)
  })
  
  if (length(cells) == 0) {
    stop("Error: No cells found for the specified cell types.")
  }
  
  # Ensure valid cell names
  valid_cells <- cells[cells %in% colnames(expression_matrix)]
  
  if (length(valid_cells) == 0) {
    stop("Error: No matching cells found in expression_matrix for the given cell types.")
  }
  
  # Check for duplicate cell names
  if (any(duplicated(valid_cells))) {
    warning("Duplicate cell names found. Using first occurrence of each cell.")
    valid_cells <- unique(valid_cells)
  }
  
  # Subset expression matrix with error handling
  expression_mat <- tryCatch({
    expression_matrix[, valid_cells, drop = FALSE]
  }, error = function(e) {
    stop("Error subsetting expression_matrix: ", e$message)
  })
  
  # Check for empty expression matrix
  if (nrow(expression_mat) == 0) {
    stop("Error: No genes found in the expression matrix.")
  }
  
  # Subset cell metadata safely
  cell_metadata_subset <- tryCatch({
    cell_metadata[match(valid_cells, rownames(cell_metadata)), , drop = FALSE]
  }, error = function(e) {
    stop("Error subsetting cell_metadata: ", e$message)
  })
  rownames(cell_metadata_subset) <- valid_cells

  # Ensure expression matrix is not empty
  if (ncol(expression_mat) == 0) {
    stop("Error: The subsetted expression matrix has zero columns. No valid cells found.")
  }

  # Ensure spatial coordinates match selected cells
  matching_spatial_cells <- intersect(valid_cells, rownames(spatial_coordinates))
  
  if (length(matching_spatial_cells) == 0) {
    stop("Error: No matching spatial coordinates found for the selected cells.")
  }

  # Subset spatial coordinates
  spatial_coordinates_subset <- tryCatch({
    spatial_coordinates[matching_spatial_cells, , drop = FALSE]
  }, error = function(e) {
    stop("Error subsetting spatial coordinates: ", e$message)
  })
  
  # Check for zero variance genes
  zero_var_genes <- apply(expression_mat, 1, function(x) var(x) == 0)
  if (any(zero_var_genes)) {
    warning("Removing ", sum(zero_var_genes), " genes with zero variance.")
    expression_mat <- expression_mat[!zero_var_genes, , drop = FALSE]
  }
  
  # Construct and return object
  object <- list(
    expression_matrix = expression_mat,
    cell_metadata = cell_metadata_subset,
    feature_metadata = feature_name,
    spatial_coordinates = data.frame(
      barcodes = matching_spatial_cells,
      x = spatial_coordinates_subset$x,
      y = spatial_coordinates_subset$y
    )
  )
  
  # Validate final object
  if (nrow(object$expression_matrix) == 0 || ncol(object$expression_matrix) == 0) {
    stop("Error: Final expression matrix is empty.")
  }
  
  if (nrow(object$cell_metadata) == 0) {
    stop("Error: Final cell metadata is empty.")
  }
  
  if (nrow(object$spatial_coordinates) == 0) {
    stop("Error: Final spatial coordinates are empty.")
  }
  
  return(object)
} 