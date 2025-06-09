#' Identify Overall Differentially Expressed Features (DEGs/DARs)
#'
#' This function identifies features that are differentially expressed in both trajectory and spatial analyses.
#'
#' @param traj_feature A character vector of trajectory-associated differentially expressed features.
#' @param sp_feature A list containing significant spatially variable features and their Moran's I values.
#'
#' @return A sorted character vector of overall significant DEGs, ordered by Moran's I deviation.
#'
#' @export
get_sig_features_geneset <- function(traj_feature, sp_feature) {
  # Input validation
  if (missing(traj_feature) || missing(sp_feature)) {
    stop("Error: Both traj_feature and sp_feature arguments are required.")
  }
  
  if (!is.list(traj_feature) || !is.list(sp_feature)) {
    stop("Error: Both traj_feature and sp_feature must be lists.")
  }
  
  required_traj_elements <- c("significant_features", "correlation_rho")
  if (!all(required_traj_elements %in% names(traj_feature))) {
    stop("Error: traj_feature must contain 'significant_features' and 'correlation_rho'.")
  }
  
  required_sp_elements <- c("significant_features")
  if (!all(required_sp_elements %in% names(sp_feature))) {
    stop("Error: sp_feature must contain 'significant_features'.")
  }
  
  if (length(traj_feature$significant_features) == 0) {
    stop("Error: No significant features found in traj_feature.")
  }
  
  if (length(sp_feature$significant_features) == 0) {
    stop("Error: No significant features found in sp_feature.")
  }
  
  # Check for missing values
  if (any(is.na(traj_feature$correlation_rho))) {
    warning("NA values found in correlation_rho. These will be removed.")
    valid_idx <- !is.na(traj_feature$correlation_rho)
    traj_feature$significant_features <- traj_feature$significant_features[valid_idx]
    traj_feature$correlation_rho <- traj_feature$correlation_rho[valid_idx]
  }

  # Find overlapping features
  overall_feature <- intersect(traj_feature$significant_features, sp_feature$significant_features)
  
  if (length(overall_feature) == 0) {
    warning("No overlapping features found between trajectory and spatial analyses.")
    return(list(
      significant_features = character(0),
      statistics = numeric(0),
      cor.rho = numeric(0)
    ))
  }

  # Get corresponding statistics
  indices_traj <- match(overall_feature, traj_feature$significant_features)
  cor_rho <- traj_feature$correlation_rho[indices_traj]

  indices_sp <- match(overall_feature, sp_feature$significant_features)
  
  # Check which statistics are available in sp_feature
  if ("Deviation_from_expectation" %in% names(sp_feature)) {
    feature_stat <- sp_feature$Deviation_from_expectation[indices_sp]
    stat_name <- "moransI"
  } else if ("importance_scores" %in% names(sp_feature)) {
    feature_stat <- sp_feature$importance_scores[indices_sp]
    stat_name <- "importance_scores"
  } else {
    feature_stat <- traj_feature$correlation_rho[indices_traj]
    stat_name <- "correlation_with_pseudotime"
  }  

  # Create and sort results
  dt <- data.frame(
    sig_feature = overall_feature,
    statistics = feature_stat,
    cor.rho = cor_rho
  )
  
  # Rename the statistics column based on the type
  names(dt)[names(dt) == "statistics"] <- stat_name
  
  dt <- dt[order(-dt[[stat_name]]), ]
  
  return(list(
    significant_features = dt$sig_feature,
    statistics = dt[[stat_name]],
    cor.rho = dt$cor.rho
  ))
}

#' Identify Overall Differentially Expressed Features (DEGs/DARs)
#'
#' This function identifies features that are differentially expressed in both trajectory and spatial analyses.
#'
#' @param obj A list containing feature metadata.
#' @param traj_feature A character vector of trajectory-associated differentially expressed features.
#' @param sp_feature A list containing significant spatially variable features and their Moran's I values.
#'
#' @return A list containing up-regulated and down-regulated features with their statistics.
#'
#' @export
get_sig_features_region <- function(obj, traj_feature, sp_feature) {
  # Input validation
  if (missing(obj) || missing(traj_feature) || missing(sp_feature)) {
    stop("Error: All arguments are required.")
  }
  
  if (!is.list(obj) || !is.list(traj_feature) || !is.list(sp_feature)) {
    stop("Error: All arguments must be lists.")
  }
  
  if (!"feature_metadata" %in% names(obj) || !"gene_short_name" %in% names(obj$feature_metadata)) {
    stop("Error: obj must contain feature_metadata with gene_short_name.")
  }
  
  required_traj_elements <- c("significant_features", "correlation_rho")
  if (!all(required_traj_elements %in% names(traj_feature))) {
    stop("Error: traj_feature must contain 'significant_features' and 'correlation_rho'.")
  }
  
  required_sp_elements <- c("significant_features", "Deviation_from_expectation")
  if (!all(required_sp_elements %in% names(sp_feature))) {
    stop("Error: sp_feature must contain 'significant_features' and 'Deviation_from_expectation'.")
  }
  
  if (length(traj_feature$significant_features) == 0) {
    stop("Error: No significant features found in traj_feature.")
  }
  
  if (length(sp_feature$significant_features) == 0) {
    stop("Error: No significant features found in sp_feature.")
  }
  
  # Check for missing values
  if (any(is.na(traj_feature$correlation_rho))) {
    warning("NA values found in correlation_rho. These will be removed.")
    valid_idx <- !is.na(traj_feature$correlation_rho)
    traj_feature$significant_features <- traj_feature$significant_features[valid_idx]
    traj_feature$correlation_rho <- traj_feature$correlation_rho[valid_idx]
  }
  
  if (any(is.na(sp_feature$Deviation_from_expectation))) {
    warning("NA values found in Deviation_from_expectation. These will be removed.")
    valid_idx <- !is.na(sp_feature$Deviation_from_expectation)
    sp_feature$significant_features <- sp_feature$significant_features[valid_idx]
    sp_feature$Deviation_from_expectation <- sp_feature$Deviation_from_expectation[valid_idx]
  }

  all_features <- obj$feature_metadata$gene_short_name
  
  if (length(all_features) == 0) {
    stop("Error: No features found in feature metadata.")
  }

  # Process trajectory features
  traj_sig_feature <- traj_feature$significant_features
  feature_corr <- traj_feature$correlation_rho
  
  if (length(traj_sig_feature) != length(feature_corr)) {
    stop("Error: Length of significant features and correlation values must match in traj_feature.")
  }
  
  dt_traj_sig_feature <- data.frame(feature = traj_sig_feature, corr = feature_corr)

  dt_traj_sig_feature_pos <- dt_traj_sig_feature[dt_traj_sig_feature$corr > 0, ]
  dt_traj_sig_feature_neg <- dt_traj_sig_feature[dt_traj_sig_feature$corr < 0, ]
  
  if (nrow(dt_traj_sig_feature_pos) == 0 && nrow(dt_traj_sig_feature_neg) == 0) {
    warning("No significant trajectory features found.")
  }
  
  dt_traj_other_features_pos <- data.frame(
    feature = all_features[!all_features %in% dt_traj_sig_feature_pos$feature],
    corr = 0
  )
  dt_traj_other_features_neg <- data.frame(
    feature = all_features[!all_features %in% dt_traj_sig_feature_neg$feature],
    corr = 0
  )

  dt_traj_pos <- rbind(dt_traj_sig_feature_pos, dt_traj_other_features_pos)
  dt_traj_neg <- rbind(dt_traj_sig_feature_neg, dt_traj_other_features_neg)

  dt_traj_pos <- dt_traj_pos[order(-abs(dt_traj_pos$corr)), ]
  dt_traj_neg <- dt_traj_neg[order(-abs(dt_traj_neg$corr)), ]

  dt_traj_pos$scaled_score_traj <- 1/seq_len(nrow(dt_traj_pos))
  dt_traj_neg$scaled_score_traj <- 1/seq_len(nrow(dt_traj_neg))

  # Process spatial features
  sp_sig_feature <- sp_feature$significant_features
  feature_moranI <- sp_feature$Deviation_from_expectation
  
  if (length(sp_sig_feature) != length(feature_moranI)) {
    stop("Error: Length of significant features and Moran's I values must match in sp_feature.")
  }
  
  dt_sp_sig_feature <- data.frame(feature = sp_sig_feature, moransI = feature_moranI)
  dt_sp_other_feature <- data.frame(
    feature = all_features[!all_features %in% sp_sig_feature],
    moransI = 0
  )
  dt_sp <- rbind(dt_sp_sig_feature, dt_sp_other_feature)
  dt_sp <- dt_sp[order(-dt_sp$moransI), ]
  dt_sp$scaled_score_sp <- 1/seq_len(nrow(dt_sp))

  # Find overlapping features
  traj_pos_features <- dt_traj_sig_feature_pos$feature
  overlap_features_pos <- intersect(traj_pos_features, sp_sig_feature)
  traj_unique_features_pos <- traj_pos_features[!traj_pos_features %in% overlap_features_pos]
  sp_unique_features_pos <- sp_sig_feature[!sp_sig_feature %in% overlap_features_pos]

  traj_neg_features <- dt_traj_sig_feature_neg$feature
  overlap_features_neg <- intersect(traj_neg_features, sp_sig_feature)
  traj_unique_features_neg <- traj_neg_features[!traj_neg_features %in% overlap_features_neg]
  sp_unique_features_neg <- sp_sig_feature[!sp_sig_feature %in% overlap_features_neg]

  # Rank aggregation with error handling
  tryCatch({
    dt_pos <- merge(dt_traj_pos, dt_sp, by = "feature")
    dt_pos <- dt_pos %>% mutate(rank_score = case_when(
      feature %in% overlap_features_pos ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp),
      feature %in% c(traj_unique_features_pos, sp_unique_features_pos) ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp),
      TRUE ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp)
    ))
    dt_pos <- dt_pos[dt_pos$feature %in% c(traj_pos_features, sp_sig_feature), ]
    dt_pos <- dt_pos[order(-dt_pos$rank_score), ]

    dt_neg <- merge(dt_traj_neg, dt_sp, by = "feature")
    dt_neg <- dt_neg %>% mutate(rank_score = case_when(
      feature %in% overlap_features_pos ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp),
      feature %in% c(traj_unique_features_neg, sp_unique_features_neg) ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp),
      TRUE ~ (abs(corr) + moransI) * (scaled_score_traj * scaled_score_sp)
    ))
    dt_neg <- dt_neg[dt_neg$feature %in% c(traj_neg_features, sp_sig_feature), ]
    dt_neg <- dt_neg[order(-dt_neg$rank_score), ]
  }, error = function(e) {
    stop("Error in rank aggregation: ", e$message)
  })

  pos_list <- list(
    significant_features = dt_pos$feature,
    statistics = dt_pos$rank_score,
    cor.rho = dt_pos$corr
  )
  neg_list <- list(
    significant_features = dt_neg$feature,
    statistics = dt_neg$rank_score,
    cor.rho = dt_neg$corr
  )

  return(list(up_regions = pos_list, down_regions = neg_list))
} 