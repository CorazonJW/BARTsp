#' Prepare Moran's I Input Data
#'
#' This function prepares the spatial and gene expression data for Moran's I calculation.
#'
#' @param object A list containing \code{expression_matrix}, \code{spatial_coordinates}, and \code{cell_metadata}.
#'
#' @return A list containing scaled expression matrix, spatial coordinates, and a formatted expression list.
#'
#' @importFrom spdep knearneigh nb2listw knn2nb moran.test
#' @importFrom stats scale p.adjust
#' @importFrom dplyr filter
#' 
#' @export
prepare_moran_input <- function(object) {
  # Input validation
  if (!is.list(object)) {
    stop("Error: object must be a list.")
  }
  
  required_elements <- c("expression_matrix", "spatial_coordinates", "cell_metadata")
  if (!all(required_elements %in% names(object))) {
    stop("Error: object must contain 'expression_matrix', 'spatial_coordinates', and 'cell_metadata'.")
  }
  
  if (!is.data.frame(object$spatial_coordinates)) {
    stop("Error: spatial_coordinates must be a data frame.")
  }
  
  if (!is.data.frame(object$cell_metadata)) {
    stop("Error: cell_metadata must be a data frame.")
  }
  
  if (!"cell_type" %in% colnames(object$cell_metadata)) {
    stop("Error: cell_metadata must contain 'cell_type' column.")
  }

  expression_matrix <- as.matrix(object$expression_matrix)
  spatial_coordinates <- object$spatial_coordinates
  cell_type_matrix <- data.frame(cell = rownames(object$cell_metadata), cell_type = object$cell_metadata$cell_type)

  features <- rownames(expression_matrix)
  if (length(features) == 0) {
    stop("Error: No features (genes) found in expression matrix.")
  }

  y_tmp <- expression_matrix
  S_tmp <- spatial_coordinates
  S_tmp <- as.matrix(S_tmp[, -1])

  # Ensure rownames(S_tmp) == colnames(y_tmp)
  if (!all(rownames(S_tmp) == colnames(y_tmp))) {
    stop("Error: Cell names in spatial coordinates do not match cell names in expression matrix.")
  }

  # Check for missing values
  if (any(is.na(y_tmp))) {
    stop("Error: Expression matrix contains missing values.")
  }
  if (any(is.na(S_tmp))) {
    stop("Error: Spatial coordinates contain missing values.")
  }

  y <- list()
  for (i in features) {
    y[[i]] <- as.matrix(y_tmp[as.character(i), ])
    colnames(y[[i]]) <- as.character(i)
  }

  return(list(expression_matrix = y_tmp, features = features, spatial_coordinates = S_tmp, y = y))
}

#' Preprocess Data for Moran's I Analysis
#'
#' This function scales gene expression and spatial coordinate matrices for Moran's I calculations.
#'
#' @param moran_object A list containing spatial data and gene expression.
#'
#' @return A modified list with centered and scaled matrices.
#'
#' @export
preprocess_data <- function(moran_object) {
  # Input validation
  if (!is.list(moran_object)) {
    stop("Error: moran_object must be a list.")
  }
  
  required_elements <- c("y", "spatial_coordinates")
  if (!all(required_elements %in% names(moran_object))) {
    stop("Error: moran_object must contain 'y' and 'spatial_coordinates'.")
  }
  
  if (!is.list(moran_object$y)) {
    stop("Error: moran_object$y must be a list.")
  }
  
  if (!is.matrix(moran_object$spatial_coordinates)) {
    stop("Error: moran_object$spatial_coordinates must be a matrix.")
  }

  scaled_y <- list()
  for (i in names(moran_object$y)) {
    # Check for non-numeric values
    if (!is.numeric(moran_object$y[[i]])) {
      stop(paste("Error: Non-numeric values found in expression data for feature:", i))
    }
    
    # Check for infinite values
    if (any(is.infinite(moran_object$y[[i]]))) {
      stop(paste("Error: Infinite values found in expression data for feature:", i))
    }
    
    scaled_y[[i]] <- scale(moran_object$y[[i]], center = TRUE, scale = TRUE)
  }

  # Check for non-numeric values in spatial coordinates
  if (!is.numeric(moran_object$spatial_coordinates)) {
    stop("Error: Non-numeric values found in spatial coordinates.")
  }
  
  # Check for infinite values in spatial coordinates
  if (any(is.infinite(moran_object$spatial_coordinates))) {
    stop("Error: Infinite values found in spatial coordinates.")
  }

  moran_object$scaled_y <- scaled_y
  moran_object$scaled_S <- scale(moran_object$spatial_coordinates, center = TRUE, scale = TRUE)
  
  return(moran_object)
}

#' Compute Moran's I for Spatial Autocorrelation
#'
#' This function computes Moran's I for each gene to assess spatial autocorrelation.
#'
#' @param moran_object A list containing scaled spatial data and gene expression.
#'
#' @return A list containing Moran's I test results.
#'
#' @importFrom spdep knearneigh nb2listw knn2nb moran.test
#' @importFrom stats scale
#'
#' @export
compute_morans_I <- function(moran_object) {
  # Input validation
  if (!is.list(moran_object)) {
    stop("Error: moran_object must be a list.")
  }
  
  required_elements <- c("scaled_S", "scaled_y")
  if (!all(required_elements %in% names(moran_object))) {
    stop("Error: moran_object must contain 'scaled_S' and 'scaled_y'.")
  }
  
  if (!is.matrix(moran_object$scaled_S)) {
    stop("Error: moran_object$scaled_S must be a matrix.")
  }
  
  if (!is.list(moran_object$scaled_y)) {
    stop("Error: moran_object$scaled_y must be a list.")
  }

  S <- moran_object$scaled_S
  coords <- as.data.frame(S)
  k <- min(5, nrow(coords) - 1)

  # Ensure there are enough points for k-nearest neighbors
  if (nrow(coords) < 2) {
    stop("Error: Not enough spatial points for nearest neighbor analysis.")
  }

  # Check for missing values
  if (any(is.na(coords))) {
    stop("Error: Missing values found in spatial coordinates.")
  }

  neighbors <- tryCatch({
    spdep::knearneigh(coords, k = k)
  }, error = function(e) {
    stop("Error computing nearest neighbors: ", e$message)
  })

  listw <- tryCatch({
    spdep::nb2listw(spdep::knn2nb(neighbors), style = "W")
  }, error = function(e) {
    stop("Error creating spatial weights: ", e$message)
  })

  scaled_y <- moran_object$scaled_y
  morans_test <- list()

  for (i in names(scaled_y)) {
    y <- scaled_y[[i]]

    # Ensure y is numeric
    if (!is.numeric(y)) {
      message(paste("Skipping feature:", i, " - y is not numeric"))
      next
    }

    # Skip if all values in y are zero
    if (all(y == 0, na.rm = TRUE)) {
      message(paste("Skipping feature:", i, " - all values are zero"))
      next
    }

    # Skip if y contains NA values
    if (any(is.na(y))) {
      message(paste("Skipping feature:", i, " - contains NA values"))
      next
    }

    # Skip if y contains infinite values
    if (any(is.infinite(y))) {
      message(paste("Skipping feature:", i, " - contains infinite values"))
      next
    }

    # Compute Moran's I
    morans_test[[i]] <- tryCatch({
      spdep::moran.test(y, listw)
    }, error = function(e) {
      message(paste("Skipping feature:", i, " - Moran's test failed:", e$message))
      return(NULL)
    })
  }

  if (length(morans_test) == 0) {
    stop("Error: No valid features for Moran's I computation.")
  }

  return(morans_test)
}

#' Extract Significant Moran's I Results
#'
#' This function extracts significant spatially autocorrelated genes from Moran's I results.
#'
#' @param morana_I_result A list of Moran's I test results.
#' @param adj.val A numeric value specifying the adjusted p-value threshold.
#' @param moransI A numeric value specifying the minimum Moran's I value.
#'
#' @return A list containing significant genes, their adjusted p-values, and deviation from expectation.
#'
#' @importFrom stats p.adjust
#'
#' @export
get_moran_result <- function(morana_I_result, adj.val, moransI) {
  # Input validation
  if (!is.list(morana_I_result)) {
    stop("Error: morana_I_result must be a list.")
  }
  
  if (length(morana_I_result) == 0) {
    stop("Error: morana_I_result is empty.")
  }
  
  if (!is.numeric(adj.val) || adj.val < 0 || adj.val > 1) {
    stop("Error: adj.val must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(moransI)) {
    stop("Error: moransI must be a numeric value.")
  }

  # Filter results
  significant_result <- tryCatch({
    morana_I_result[sapply(morana_I_result, function(x) {
      !is.null(x) && x$estimate[[1]] != x$estimate[[2]]
    })]
  }, error = function(e) {
    stop("Error filtering results: ", e$message)
  })

  significant_result <- significant_result[sapply(significant_result, function(x) {
    x$adjusted_p.value < as.numeric(adj.val)
  })]

  significant_result <- significant_result[sapply(significant_result, function(x) {
    x$estimate[[1]] >= x$estimate[[2]] + moransI
  })]

  if (length(significant_result) == 0) {
    warning("No significant results found at the specified thresholds.")
    return(list(
      significant_features = character(0),
      adjusted_pvals = numeric(0),
      Deviation_from_expectation = numeric(0)
    ))
  }

  adjusted_pvals <- sapply(significant_result, function(x) x$adjusted_p.value)
  names(adjusted_pvals) <- NULL
  deviation <- sapply(significant_result, function(x) x$estimate[[1]] - x$estimate[[2]])
  names(deviation) <- NULL

  return(list(
    significant_features = names(significant_result),
    adjusted_pvals = adjusted_pvals,
    Deviation_from_expectation = deviation
  ))
}