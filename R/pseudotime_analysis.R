#' Construct a Trajectory Using Monocle3
#'
#' @param object A list containing expression matrix, cell metadata, and gene metadata.
#' @param start_cell_type A character string indicating the starting cell type for trajectory inference.
#'
#' @return A \code{monocle3} cell_data_set object with trajectory information.
#'
#' @importFrom monocle3 new_cell_data_set preprocess_cds reduce_dimension cluster_cells learn_graph order_cells
#' @importFrom SeuratWrappers some_function_if_used
#' @importFrom igraph V
#' 
#' @export
construct_trajectory <- function(object, start_cell_type) {
  # Input validation
  if (!is.list(object)) {
    stop("Error: object must be a list.")
  }
  
  required_elements <- c("expression_matrix", "cell_metadata", "feature_metadata")
  if (!all(required_elements %in% names(object))) {
    stop("Error: object must contain 'expression_matrix', 'cell_metadata', and 'feature_metadata'.")
  }
  
  if (!is.character(start_cell_type) || length(start_cell_type) != 1) {
    stop("Error: start_cell_type must be a single character string.")
  }
  
  if (!start_cell_type %in% object$cell_metadata$cell_type) {
    stop("Error: start_cell_type not found in cell_metadata$cell_type.")
  }
  
  # Check for missing values
  if (any(is.na(object$expression_matrix))) {
    stop("Error: expression_matrix contains missing values.")
  }
  
  if (any(is.na(object$cell_metadata$cell_type))) {
    stop("Error: cell_metadata contains missing values in cell_type column.")
  }
  
  # Check for zero variance genes
  zero_var_genes <- apply(object$expression_matrix, 1, function(x) var(x) == 0)
  if (any(zero_var_genes)) {
    warning("Removing ", sum(zero_var_genes), " genes with zero variance.")
    object$expression_matrix <- object$expression_matrix[!zero_var_genes, , drop = FALSE]
  }

  expression_matrix <- object$expression_matrix
  cell_metadata <- object$cell_metadata
  feature_metadata <- object$feature_metadata

  # Create monocle3 object with error handling
  cds <- tryCatch({
    monocle3::new_cell_data_set(expression_matrix, 
                               cell_metadata = cell_metadata, 
                               feature_metadata = feature_metadata)
  }, error = function(e) {
    stop("Error creating cell_data_set: ", e$message)
  })

  # Preprocess with error handling
  cds <- tryCatch({
    monocle3::preprocess_cds(cds, num_dim = 50)
  }, error = function(e) {
    stop("Error in preprocessing: ", e$message)
  })
  
  cds <- tryCatch({
    monocle3::reduce_dimension(cds)
  }, error = function(e) {
    stop("Error in dimension reduction: ", e$message)
  })

  # Learn trajectory graph with error handling
  cds <- tryCatch({
    monocle3::cluster_cells(cds)
  }, error = function(e) {
    stop("Error in cell clustering: ", e$message)
  })
  
  cds <- tryCatch({
    monocle3::learn_graph(cds)
  }, error = function(e) {
    stop("Error in graph learning: ", e$message)
  })

  # Order cells based on pseudotime
  get_earliest_principal_node <- function(cds, time_bin = start_cell_type) {
    cell_ids <- which(SummarizedExperiment::colData(cds)[, "cell_type"] == time_bin)
    
    if (length(cell_ids) == 0) {
      stop("Error: No cells found for the specified start_cell_type.")
    }

    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    if (is.null(closest_vertex)) {
      stop("Error: Principal graph not found. Check if dimension reduction was successful.")
    }
    
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    
    root_pr_nodes <- tryCatch({
      igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[
        as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
      ]
    }, error = function(e) {
      stop("Error finding root nodes: ", e$message)
    })
    
    if (length(root_pr_nodes) == 0) {
      stop("Error: Could not determine root nodes for trajectory.")
    }
    
    root_pr_nodes
  }

  cds <- tryCatch({
    monocle3::order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
  }, error = function(e) {
    stop("Error ordering cells: ", e$message)
  })

  return(cds)
}

#' Identify Trajectory-Associated Differentially Expressed Features (DEGs/DARs)
#'
#' This function identifies genes whose expression is correlated with pseudotime along a trajectory.
#'
#' @param pseudotime A numeric vector of pseudotime values.
#' @param object A list containing \code{expression_matrix} (expression data) 
#' and \code{cell_metadata} (cell annotations).
#' @param pval_cutoff A numeric value specifying the adjusted p-value threshold for significance.
#' @param cor_cutoff_pos A numeric value specifying the positive correlation threshold.
#' @param cor_cutoff_neg A numeric value specifying the negative correlation threshold.
#'
#' @return A list containing significant features, adjusted p-values, and correlation values.
#'
#' @importFrom stats cor.test p.adjust
#'
#' @export
get_traj_features <- function(pseudotime, object, pval_cutoff, cor_cutoff_pos, cor_cutoff_neg) {
  # Input validation
  if (missing(pseudotime) || missing(object) || missing(pval_cutoff) || 
      missing(cor_cutoff_pos) || missing(cor_cutoff_neg)) {
    stop("Error: All arguments are required.")
  }
  
  if (!is.numeric(pseudotime)) {
    stop("Error: pseudotime must be a numeric vector.")
  }
  
  if (!is.list(object)) {
    stop("Error: object must be a list.")
  }
  
  if (!"expression_matrix" %in% names(object)) {
    stop("Error: object must contain 'expression_matrix'.")
  }
  
  if (!is.numeric(pval_cutoff) || pval_cutoff < 0 || pval_cutoff > 1) {
    stop("Error: pval_cutoff must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(cor_cutoff_pos) || cor_cutoff_pos < 0 || cor_cutoff_pos > 1) {
    stop("Error: cor_cutoff_pos must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(cor_cutoff_neg) || cor_cutoff_neg < -1 || cor_cutoff_neg > 0) {
    stop("Error: cor_cutoff_neg must be a numeric value between -1 and 0.")
  }

  # Extract data
  expression_mat <- object$expression_matrix
  
  # Check for missing values
  if (any(is.na(expression_mat))) {
    stop("Error: expression_matrix contains missing values.")
  }
  
  if (any(is.na(pseudotime))) {
    warning("Pseudotime contains missing values. These will be removed.")
  }

  # Ensure column names and pseudotime names are character
  colnames(expression_mat) <- as.character(colnames(expression_mat))
  names(pseudotime) <- as.character(names(pseudotime))

  # Extract pseudotime values
  pseudotime_values <- setNames(as.numeric(pseudotime), names(pseudotime))

  # Ensure pseudotime matches expression matrix
  common_cells <- intersect(colnames(expression_mat), names(pseudotime_values))

  if (length(common_cells) == 0) {
    stop("Error: No matching cells between expression matrix and pseudotime data.")
  }

  # Subset data to common cells
  expression_mat <- expression_mat[, common_cells, drop=FALSE]
  pseudotime_values <- pseudotime_values[common_cells, drop=FALSE]

  # Ensure expression matrix is not empty after subsetting
  if (ncol(expression_mat) == 0) {
    stop("Error: Expression matrix has zero columns after filtering with common cells.")
  }

  # Remove NAs from pseudotime
  valid_cells <- !is.na(pseudotime_values)
  expression_mat <- expression_mat[, valid_cells, drop = FALSE]
  pseudotime_values <- pseudotime_values[valid_cells]

  if (length(pseudotime_values) == 0) {
    stop("Error: All pseudotime values are NA.")
  }

  # Check for zero variance genes
  zero_var_genes <- apply(expression_mat, 1, function(x) var(x) == 0)
  if (any(zero_var_genes)) {
    warning("Removing ", sum(zero_var_genes), " genes with zero variance.")
    expression_mat <- expression_mat[!zero_var_genes, , drop = FALSE]
  }

  # Correlate each gene with pseudotime
  cor_results <- tryCatch({
    apply(expression_mat, 1, function(gene_expr) {
      if (all(is.na(gene_expr))) return(NA)  # Avoid errors on empty genes
      test <- stats::cor.test(gene_expr, pseudotime_values, method = "spearman")
      c(correlation = test$estimate, p_value = test$p.value)
    })
  }, error = function(e) {
    stop("Error computing correlations: ", e$message)
  })
  
  # Adjust p-values using FDR correction
  cor_results <- t(cor_results) 
  colnames(cor_results) <- c("correlation.rho", "p_value")
  cor_results_df <- as.data.frame(cor_results)
  
  # Check for NA values in results
  if (any(is.na(cor_results_df$p_value))) {
    warning("Some genes had NA p-values and will be removed.")
  }
  
  cor_results_df$adjusted_pval <- stats::p.adjust(cor_results_df$p_value, method = "fdr")
  cor_results_df <- cor_results_df[!is.na(cor_results_df$adjusted_pval), ]

  # Extract significant genes
  dt <- cor_results_df[cor_results_df$adjusted_pval < pval_cutoff &
                     ((cor_results_df$correlation.rho > cor_cutoff_pos & cor_results_df$correlation.rho > 0) |
                      (cor_results_df$correlation.rho < cor_cutoff_neg & cor_results_df$correlation.rho < 0)), ]

  if (nrow(dt) == 0) {
    warning("No significant features found at the specified thresholds.")
    return(list(
      significant_features = character(0),
      adjusted_pvals = numeric(0),
      correlation_rho = numeric(0)
    ))
  }

  dt <- dt[order(-abs(dt$correlation.rho)), ]

  significant_features <- rownames(dt) 
  adjusted_pvals <- dt$adjusted_pval
  correlation_rho <- dt$correlation.rho

  return(list(
    significant_features = significant_features,
    adjusted_pvals = adjusted_pvals,
    correlation_rho = correlation_rho
  ))
} 