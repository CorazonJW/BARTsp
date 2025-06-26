#' Integrate BART Upstream and Downstream Results into a Combined Rank
#'
#' This function integrates upstream and downstream BART results, filters TFs by statistical significance, 
#' and computes a combined rank score based on scaled ranks. TFs with either upstream or downstream 
#' significance (based on user-defined p-value cutoffs) are retained and re-ranked.
#'
#' @param results_geneset_up Data frame containing upstream TF results with columns `TF`, `rank_avg_z_p_a_irwinhall_pvalue`.
#' @param results_geneset_down Data frame containing downstream TF results with the same structure as upstream.
#' @param cutoff_up Numeric. P-value threshold for upstream significance.
#' @param cutoff_down Numeric. P-value threshold for downstream significance.
#'
#' @return A data frame with TF names, new rank order, and combined rank score.
#' @export
integrate_bart_result <- function(results_geneset_up, results_geneset_down, cutoff_up, cutoff_down) {
  # Input validation
  if (missing(results_geneset_up) || missing(results_geneset_down) || missing(cutoff_up) || missing(cutoff_down)) {
    stop("Error: All arguments (results_geneset_up, results_geneset_down, cutoff_up, cutoff_down) are required.")
  }

  if (!is.data.frame(results_geneset_up) || !is.data.frame(results_geneset_down)) {
    stop("Error: Input results must be data frames.")
  }

  if (!is.numeric(cutoff_up) || !is.numeric(cutoff_down)) {
    stop("Error: Cutoff values must be numeric.")
  }

  if (!all(c("TF", "rank_avg_z_p_a_irwinhall_pvalue") %in% colnames(results_geneset_up))) {
    stop("Error: results_geneset_up must contain columns 'TF' and 'rank_avg_z_p_a_irwinhall_pvalue'.")
  }

  if (!all(c("TF", "rank_avg_z_p_a_irwinhall_pvalue") %in% colnames(results_geneset_down))) {
    stop("Error: results_geneset_down must contain columns 'TF' and 'rank_avg_z_p_a_irwinhall_pvalue'.")
  }

  # Assign ranks if not already present
  results_geneset_up$rank <- as.numeric(seq_len(nrow(results_geneset_up)))
  results_geneset_down$rank <- as.numeric(seq_len(nrow(results_geneset_down)))

  # Extract TF names and check overlap
  upstream_TF <- results_geneset_up$TF
  downstream_TF <- results_geneset_down$TF

  if (!setequal(upstream_TF, downstream_TF)) {
    stop("Error: TFs in upstream and downstream results do not match.")
  }

  # Construct rank and p-value data frames
  upstream_TF_dt <- data.frame(TF = results_geneset_up$TF,
                               upstream_rank = results_geneset_up$rank,
                               upstream_pvalue = results_geneset_up$rank_avg_z_p_a_irwinhall_pvalue,
                               upstream_scaled_rank = 1 / results_geneset_up$rank)

  downstream_TF_dt <- data.frame(TF = results_geneset_down$TF,
                                 downstream_rank = results_geneset_down$rank,
                                 downstream_pvalue = results_geneset_down$rank_avg_z_p_a_irwinhall_pvalue,
                                 downstream_scaled_rank = 1 / results_geneset_up$rank)

  # Merge and calculate combined rank score
  dt <- merge(upstream_TF_dt, downstream_TF_dt, by = "TF")

  dt$TF_type <- ifelse(dt$upstream_pvalue < cutoff_up | dt$downstream_pvalue < cutoff_down, "either_sig", "not_sig")

  dt$new_rank_score <- ifelse(dt$TF_type == "either_sig", dt$upstream_scaled_rank - dt$downstream_scaled_rank, NA)

  dt <- dt[!is.na(dt$new_rank_score), ]
  dt <- dt[order(-dt$new_rank_score), ]
  dt$new_rank <- seq_len(nrow(dt))

  return(data.frame(TF = dt$TF, new_rank = dt$new_rank, rank_score = dt$new_rank_score))
}


#' Bar Plot of Top and Bottom TFs by Rank Score
#'
#' This function generates a bar plot of the top and bottom transcription factors (TFs) based on rank score 
#' from the aggregated BART results. Positive and negative rank scores are visualized in different colors.
#'
#' @param dt A data frame containing at least columns `TF` and `rank_score`.
#' @param top_n Integer. Number of top and bottom TFs to include (default = 10).
#'
#' @return A `ggplot2` object showing the ranked bar plot of TFs.
#' @export
plot_integration_bar <- function(dt, top_n = 10) {

  if (!is.data.frame(dt) || !all(c("TF", "rank_score") %in% colnames(dt))) {
    stop("Error: Input must be a data frame with columns 'TF' and 'rank_score'.")
  }

  if (!is.numeric(top_n) || top_n <= 0) {
    stop("Error: top_n must be a positive integer.")
  }

  dt_top <- dt %>% arrange(desc(rank_score)) %>% slice(1:top_n)
  dt_bottom <- dt %>% arrange(rank_score) %>% slice(1:top_n)

  dt_plot <- bind_rows(dt_top, dt_bottom) %>%
    mutate(color = ifelse(rank_score >= 0, "positive", "negative"))

  ggplot(dt_plot, aes(x = reorder(TF, rank_score), y = rank_score, fill = color)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("positive" = "#FA7E7A", "negative" = "#4EA72E")) +
    ylim(-1, 1) +
    theme_bw() +
    geom_vline(xintercept = top_n + 0.5, linetype = "dotted", color = "grey", linewidth = 1) +
    labs(x = NULL, y = "Rank Score") +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}


#' Dot Plot of All TF Rank Scores with Highlighted TFs
#'
#' Generates a dot plot for all TFs ranked by their score, optionally highlighting a user-defined subset of TFs. 
#' Highlighted TFs are annotated with labels and colored distinctly.
#'
#' @param dt A data frame containing at least columns `TF` and `rank_score`.
#' @param tf_highlight Character vector. TFs to be highlighted and labeled on the plot.
#'
#' @return A `ggplot2` object showing the dot plot with highlighted TFs.
#' @export
plot_integration_dot <- function(dt, tf_highlight = tf_interest) {

  if (!is.data.frame(dt) || !all(c("TF", "rank_score") %in% colnames(dt))) {
    stop("Error: Input must be a data frame with columns 'TF' and 'rank_score'.")
  }

  if (!is.character(tf_highlight)) {
    stop("Error: tf_highlight must be a character vector of TF names.")
  }

  dt <- dt %>% arrange(rank_score) %>% distinct(TF, .keep_all = TRUE)
  dt$TF_ordered <- factor(dt$TF, levels = rev(dt$TF))
  dt$highlight <- dt$TF %in% tf_highlight
  dt$color <- ifelse(dt$highlight, "highlight", ifelse(dt$rank_score >= 0, "positive", "negative"))

  ggplot() +
    geom_point(data = dt, aes(x = TF_ordered, y = rank_score, color = color), size = 2) +
    geom_point(data = dt %>% filter(highlight), aes(x = TF_ordered, y = rank_score), color = "#E41A1C", size = 2) +
    ggrepel::geom_label_repel(data = dt %>% filter(highlight),
                              aes(x = TF_ordered, y = rank_score, label = TF),
                              size = 3, color = "black", fill = "white", box.padding = 0.3, max.overlaps = Inf) +
    scale_color_manual(values = c("positive" = "#FA7E7A", "negative" = "#4EA72E", "highlight" = "#E41A1C")) +
    ylim(-1, 1) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    labs(x = NULL, y = "Rank Score") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}
