#' Compute Neighbors for Spatial Analysis
#'
#' This function computes the k-nearest neighbors for each cell based on spatial coordinates.
#'
#' @param spatial_coordinates A data frame containing spatial coordinates (x, y) for each cell.
#' @param k Number of nearest neighbors to compute. Default is 5.
#'
#' @return A list containing:
#' \item{neighbors}{A matrix where each row represents a cell and columns contain indices of its k-nearest neighbors}
#' \item{distances}{A matrix containing the distances to the k-nearest neighbors}
#'
#' @importFrom stats dist
#' @export
compute_neighbors <- function(spatial_coordinates, k = 5) {
  # Ensure input is valid
  if (!is.data.frame(spatial_coordinates)) {
    stop("Error: spatial_coordinates must be a data frame.")
  }
  
  if (!all(c("x", "y") %in% colnames(spatial_coordinates))) {
    stop("Error: spatial_coordinates must contain 'x' and 'y' columns.")
  }
  
  if (!is.numeric(k) || k < 1) {
    stop("Error: k must be a positive integer.")
  }
  
  # Check for missing values
  if (any(is.na(spatial_coordinates$x)) || any(is.na(spatial_coordinates$y))) {
    stop("Error: Spatial coordinates contain missing values.")
  }
  
  # Check for non-numeric values
  if (!is.numeric(spatial_coordinates$x) || !is.numeric(spatial_coordinates$y)) {
    stop("Error: Spatial coordinates must be numeric.")
  }
  
  # Check for infinite values
  if (any(is.infinite(spatial_coordinates$x)) || any(is.infinite(spatial_coordinates$y))) {
    stop("Error: Spatial coordinates contain infinite values.")
  }

  # Calculate pairwise distances
  dist_matrix <- tryCatch({
    as.matrix(stats::dist(spatial_coordinates[, c("x", "y")]))
  }, error = function(e) {
    stop("Error computing distances: ", e$message)
  })
  
  # Find k-nearest neighbors for each cell
  k <- min(k, nrow(spatial_coordinates) - 1)  # Ensure k is not larger than n-1
  neighbors <- matrix(0, nrow = nrow(spatial_coordinates), ncol = k)
  distances <- matrix(0, nrow = nrow(spatial_coordinates), ncol = k)
  
  for (i in seq_len(nrow(spatial_coordinates))) {
    # Get distances for current cell
    cell_distances <- dist_matrix[i, ]
    # Find k-nearest neighbors (excluding self)
    sorted_indices <- order(cell_distances)[2:(k + 1)]  # Skip first (self)
    neighbors[i, ] <- sorted_indices
    distances[i, ] <- cell_distances[sorted_indices]
  }
  
  return(list(neighbors = neighbors, distances = distances))
}

#' Calculate Importance Score for Spatial Genes
#'
#' This function calculates an importance score for each gene based on its spatial expression pattern.
#'
#' @param expression_matrix A matrix of gene expression values.
#' @param neighbors A matrix of k-nearest neighbors for each cell.
#' @param distances A matrix of distances to k-nearest neighbors.
#' @param method Method to calculate importance score. Options: "correlation" (default) or "variance".
#'
#' @return A data frame containing gene names and their importance scores.
#'
#' @importFrom stats cor
#' @export
calculate_importance_score <- function(expression_matrix, neighbors, distances, method = "correlation") {
  # Ensure inputs are valid
  if (!is.matrix(expression_matrix)) {
    stop("Error: expression_matrix must be a matrix.")
  }
  
  if (!is.matrix(neighbors) || !is.matrix(distances)) {
    stop("Error: neighbors and distances must be matrices.")
  }
  
  if (ncol(expression_matrix) != nrow(neighbors)) {
    stop("Error: Number of cells in expression_matrix must match number of rows in neighbors matrix.")
  }
  
  if (!method %in% c("correlation", "variance")) {
    stop("Error: method must be either 'correlation' or 'variance'.")
  }
  
  # Check for missing values
  if (any(is.na(expression_matrix))) {
    stop("Error: Expression matrix contains missing values.")
  }
  
  if (any(is.na(neighbors)) || any(is.na(distances))) {
    stop("Error: Neighbors or distances matrices contain missing values.")
  }
  
  # Check for non-numeric values
  if (!is.numeric(expression_matrix)) {
    stop("Error: Expression matrix must contain numeric values.")
  }
  
  if (!is.numeric(neighbors) || !is.numeric(distances)) {
    stop("Error: Neighbors and distances matrices must contain numeric values.")
  }
  
  # Check for infinite values
  if (any(is.infinite(expression_matrix))) {
    stop("Error: Expression matrix contains infinite values.")
  }
  
  if (any(is.infinite(neighbors)) || any(is.infinite(distances))) {
    stop("Error: Neighbors or distances matrices contain infinite values.")
  }

  # Initialize results
  n_genes <- nrow(expression_matrix)
  importance_scores <- numeric(n_genes)
  
  # Calculate importance score for each gene
  for (i in seq_len(n_genes)) {
    gene_expr <- expression_matrix[i, ]
    
    if (method == "correlation") {
      # Calculate correlation between gene expression and spatial distances
      spatial_scores <- numeric(ncol(neighbors))
      for (j in seq_len(ncol(neighbors))) {
        neighbor_expr <- gene_expr[neighbors[, j]]
        spatial_scores[j] <- mean(abs(stats::cor(gene_expr, neighbor_expr, method = "spearman")))
      }
      importance_scores[i] <- mean(spatial_scores)
    } else if (method == "variance") {
      # Calculate variance of expression in spatial neighborhoods
      spatial_variances <- numeric(ncol(neighbors))
      for (j in seq_len(ncol(neighbors))) {
        neighbor_expr <- gene_expr[neighbors[, j]]
        spatial_variances[j] <- var(neighbor_expr)
      }
      importance_scores[i] <- mean(spatial_variances)
    }
  }
  
  # Create results data frame
  results <- data.frame(
    gene = rownames(expression_matrix),
    importance_score = importance_scores,
    stringsAsFactors = FALSE
  )
  
  # Sort by importance score
  results <- results[order(-results$importance_score), ]
  
  return(results)
}

#' Run KNN Spatial Analysis
#'
#' This function performs spatial gene analysis using k-nearest neighbors and importance scores.
#'
#' @param object A list containing \code{expression_matrix} and \code{spatial_coordinates}.
#' @param k Number of nearest neighbors to consider. Default is 5.
#' @param method Method to calculate importance score. Options: "correlation" (default) or "variance".
#' @param cutoff Importance score threshold for significance.
#'
#' @return A list containing:
#' \item{significant_genes}{Character vector of significant spatially variable genes}
#' \item{importance_scores}{Data frame of all genes and their importance scores}
#' \item{neighbors}{Matrix of k-nearest neighbors}
#'
#' @export
run_knn_spatial <- function(object, k = 5, method = "correlation", cutoff = NULL) {
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
  
  if (!is.null(cutoff) && (!is.numeric(cutoff) || cutoff < 0)) {
    stop("Error: cutoff must be a non-negative numeric value.")
  }
  
  # Check for missing values
  if (any(is.na(object$expression_matrix))) {
    stop("Error: Expression matrix contains missing values.")
  }
  
  if (any(is.na(object$spatial_coordinates))) {
    stop("Error: Spatial coordinates contain missing values.")
  }
  
  # Compute neighbors
  neighbors_result <- tryCatch({
    compute_neighbors(object$spatial_coordinates, k = k)
  }, error = function(e) {
    stop("Error computing neighbors: ", e$message)
  })
  
  # Calculate importance scores
  importance_scores <- tryCatch({
    calculate_importance_score(
      object$expression_matrix,
      neighbors_result$neighbors,
      neighbors_result$distances,
      method = method
    )
  }, error = function(e) {
    stop("Error calculating importance scores: ", e$message)
  })
  
  # Get significant genes if cutoff is provided
  significant_genes <- NULL
  if (!is.null(cutoff)) {
    significant_genes <- importance_scores$gene[importance_scores$importance_score >= cutoff]
    if (length(significant_genes) == 0) {
      warning("No significant genes found at the specified cutoff.")
    }
  }
  
  return(list(
    significant_genes = significant_genes,
    importance_scores = importance_scores$importance_score
  ))
} 