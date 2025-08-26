
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ArchR)
library(grid)

source("~/manuscript_figure.r")

##### 1. Load SVFs detected by each method
BART_DAR <- read.csv("~/spatially_significant_peaks_BART.txt") %>% filter(adjusted_pvals < 0.05 & Deviation_from_expectation > 0.1)
BART_peak <- BART_DAR$significant_features #581
MERINGUE_result <- read.table("~/spatially_significant_peaks_MERINGUE.txt")
MERINGUE_peak <- rownames(MERINGUE_result) #158
SpaGene_result <- read.table("~/spatially_significant_peaks_SpaGene.txt")
SpaGene_peak <- rownames(SpaGene_result) #38
SPADE_result <- read.table("~/spatially_significant_peaks_SPADE.txt") %>% filter(Adjust.Pvalue < 0.1)
SPADE_peak <- SPADE_result$geneid #1752, 0 after p-adj (not very well) 
haystack_result <- read.table("~/spatially_significant_peaks_haystack.txt")
haystack_peak <- rownames(haystack_result) #1109, 0 after p-adj (not very well)
spatrack_result <- read.table("~/spatially_significant_peaks_spatrack.tsv", sep="\t", header = TRUE)
spatrack_peak_down <- spatrack_result[spatrack_result$pattern == "increase", ]$gene #30
spatrack_peak_up <- spatrack_result[spatrack_result$pattern == "decrease", ]$gene #39
spatrack_peak <- c(spatrack_peak_up, spatrack_peak_down)





##### 2. Check number of SVFs detected by each method
# Figure 6A
dt <- data.frame(method = c("BART-spatial", "MERINGUE", "SpaTrack", "SpaGene", "SPADE", "singleCell\nHaystack"), 
                 number_of_SVPs = c(length(BART_peak), length(MERINGUE_peak), length(spatrack_peak), length(SpaGene_peak), length(SPADE_peak), length(haystack_peak)))

dt <- dt[order(dt$number_of_SVPs), ]
dt$method <- factor(dt$method, levels = dt$method[order(dt$number_of_SVPs)])

set.seed(1)
palette <- c(sci_colors(4, effect = "contrast"), "#4C78A8", "#F58518")
names(palette) <- c("MERINGUE", "SpaGene", "SPADE", "singleCell\nHaystack", "BART-spatial", "SpaTrack")

p <- ggplot(dt, aes(x = number_of_SVPs, y = method, fill = method)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = number_of_SVPs), vjust = 0.5, size = 4) +
      xlim(0, 620) +
      theme_bw() +
      scale_fill_manual(values = palette) +
      labs(x = "", y = "", title = "Number of SVFs") +
      theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, family = "Arial"),
            axis.text.y = element_text(size = 12, family = "Arial"), 
            axis.title = element_blank(), 
            plot.title = element_text(size = 14, hjust = 0.5, family = "Arial", face = "bold"),
            legend.position = "none")
ggsave("svp_barplot.png", plot = p, width = 3.3, height = 3)

# Figure S3D
moran_DAR <- read.csv("~/spATAC/sp_regionset_all_peaks.csv") %>% filter(adjusted_pvals < 0.05)
temporal_DAR <- read.csv("~/tmp_regionset.csv")
common <- intersect(moran_DAR$significant_features, temporal_DAR$significant_features)
bart_up <- temporal_DAR %>% filter(significant_features %in% common & correlation_rho > 0)
bart_down <- temporal_DAR %>% filter(significant_features %in% common & correlation_rho < 0)

spatrack_result <- read.table("~/ATAC_spatrack_sort_exp_sig_with_trend.tsv", sep="\t", header = TRUE)
spatrack_up <- spatrack_result %>% filter(pattern == "increase")
spatrack_down <- spatrack_result %>% filter(pattern == "decrease")

n_up <- n_distinct(spatrack_up$gene)
n_down <- n_distinct(spatrack_down$gene)
n_up_bart <- n_distinct(bart_up$significant_features)
n_down_bart <- n_distinct(bart_down$significant_features)

tf_counts <- data.frame(Method = c("BART-spatial", "SpaTrack", "BART-spatial", "SpaTrack"), Pattern = c("decrease", "decrease", "increase", "increase"),
                        Count = c(n_down_bart, n_down, n_up_bart, n_up))

tf_counts <- tf_counts %>% mutate(x_axis = paste(Method, Pattern, sep = "_")) %>%
                           mutate(x_axis = factor(x_axis, levels = c("BART-spatial_decrease", "SpaTrack_decrease", "BART-spatial_increase", "SpaTrack_increase")))

p <- ggplot(tf_counts, aes(x = x_axis, y = Count, fill = Pattern)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, family = "Arial") +
  theme_bw() +
  ylim(0, 75) +
  labs(x = NULL, y = "Number of SVFs") +
  scale_fill_manual(values = c("decrease" = "#E3B251", "increase" = "#994C00")) +
  scale_x_discrete(labels = c("BART-\nspatial", "SpaTrack", "BART-\nspatial", "SpaTrack")) +
  theme(legend.position = "top", 
        axis.text = element_text(size = 12, family = "Arial"), 
        axis.title = element_text(size = 14, face = "bold", family = "Arial"), 
        legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
        legend.text = element_text(size = 12, family = "Arial"))

ggsave("number_of_SVFs_with_trend.png", p, height = 3, width = 4)





##### 3. Check average accessibility of SVFs detected by each method
subset <- readRDS("~/spRNA_subset.rds")

expression_matrix <- subset@assays$peaks@counts
cell_metadata <- subset@meta.data
spatial_coordinates <- data.frame(subset@images$slice1@coordinates)
spatial_coordinates$x <- spatial_coordinates$imagerow
spatial_coordinates$y <- spatial_coordinates$imagecol

dt_1 <- intersect(BART_peak, MERINGUE_peak)
dt_2 <- intersect(SpaGene_peak, BART_peak)
dt_3 <- intersect(spatrack_peak, BART_peak)
unique_bart_peak <- BART_peak[!BART_peak %in% c(dt_1, dt_2, dt_3)]

cell_df <- cell_metadata
cell_df$x <- spatial_coordinates$x 
cell_df$y <- spatial_coordinates$y

p <- peak_avg_acc(expression_matrix, unique_bart_peak, cell_df, "BART-spatial unique SVFs\naverage accessibility")
p1 <- peak_avg_acc(expression_matrix, BART_peak, cell_df, "BART-spatial SVFs\naverage accessibility")
p2 <- peak_avg_acc(expression_matrix, MERINGUE_peak, cell_df, "MERINGUE SVFs\naverage accessibility")
p3 <- peak_avg_acc(expression_matrix, spatrack_peak_up, cell_df, "SpaTrack (increase) SVFs\naverage accessibility")
p4 <- peak_avg_acc(expression_matrix, spatrack_peak_down, cell_df, "SpaTrack (decrease) SVFs\naverage accessibility")
p5 <- peak_avg_acc(expression_matrix, SpaGene_peak, cell_df, "SpaGene SVFs\naverage accessibility")
p6 <- peak_avg_acc(expression_matrix, haystack_peak, cell_df, "singleCellHaystack SVFs\naverage accessibility")
ggsave("accessibility_plots_BARTspatial_unique.png", p, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_BARTspatial.png", p1, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_MERINGUE.png", p2, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_SpaTrack_increase.png", p3, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_SpaTrack_decrease.png", p4, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_SpaGene.png", p5, height = 3, width = 4, dpi = 600)
ggsave("accessibility_plots_singleCellHaystack.png", p6, height = 3, width = 4, dpi = 600)





##### 4. Check difference of accessibility between cell types 
subset <- readRDS("/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/E13/spRNA_output/spRNA_subset.rds")

BARTspatial <- SVF_avg_acc("BART-\nspatial", cell_df, BART_peak)
MERINGUE <- SVF_avg_acc("MERINGUE\n", cell_df, MERINGUE_peak)
spatrack_up <- SVF_avg_acc("SpaTrack\n(inrease)", cell_df, spatrack_peak_up)
spatrack_down <- SVF_avg_acc("SpaTrack\n(decrease)", cell_df, spatrack_peak_down)
SpaGene <- SVF_avg_acc("SpaGene\n", cell_df, SpaGene_peak)
haystack <- SVF_avg_acc("Singlecell\nHaystack", cell_df, haystack_peak)
BARTspatial_u <- SVF_avg_acc("BART-spatial\nunique", cell_df, unique_bart_peak)

combined_df <- bind_rows(
  BARTspatial %>% select(avg_acc, method, cell_type),
  MERINGUE %>% select(avg_acc, method, cell_type), 
  spatrack_down %>% select(avg_acc, method, cell_type),
  SpaGene %>% select(avg_acc, method, cell_type)
)

p7 <- SVF_avg_acc_plot(combined_df, "two.sided")
ggsave("difference_in_SVF_acc.png", p7, height = 3, width = 7.5, dpi = 600)
p8 <- SVF_avg_acc_plot(BARTspatial_u, "two.sided")
ggsave("difference_in_SVF_acc_bart_unique.png", p8, height = 3, width = 5, dpi = 600)





##### 5. Check correlation of SVFs accessibility between cell types and pseudotime
## Compute average accessibility of SVFs
get_avg_acc <- function(method, cell_df, data) {
  moran_expr <- Matrix::colMeans(expression_matrix[data, , drop = FALSE])
  cell_df$avg_acc <- moran_expr
  cell_df$method <- as.character(method)
  return(cell_df)
}

BARTspatial <- get_avg_acc("BART-spatial", cell_df, BART_peak)
MERINGUE <- get_avg_acc("MERINGUE", cell_df, MERINGUE_peak)
spatrack_up <- get_avg_acc("SpaTrack (inrease)", cell_df, spatrack_peak_up)
spatrack_down <- get_avg_acc("SpaTrack (decrease)", cell_df, spatrack_peak_down)
SpaGene <- get_avg_acc("SpaGene", cell_df, SpaGene_peak)
haystack <- get_avg_acc("SinglecellHaystack", cell_df, haystack_peak)
BARTspatial_u <- get_avg_acc("BART-spatial unique", cell_df, unique_bart_peak)


## Steiger's Z-test
# 0) Variables
ptime <- subset$spATAC_traj
base <- list("BART-spatial" = BARTspatial, "MERINGUE" = MERINGUE, "SpaGene" = SpaGene,
             "SpaTrack (decrease)" = spatrack_down, "singleCellHaystack" = haystack, "BART-spatial unique" = BARTspatial_u)
method <- list("BART-spatial" = BARTspatial, "MERINGUE" = MERINGUE, "SpaGene" = SpaGene,
               "SpaTrack (decrease)" = spatrack_down, "singleCellHaystack" = haystack, "BART-spatial unique" = BARTspatial_u)

# 1) One-sided Steiger (pooled Z*) per PASS
steiger_one_sided <- function(base, method, n = 246, alternative = c("method>baseline","baseline>method")) {
  alternative <- match.arg(alternative)

  # correlations (pairwise handling w/o explicit mask)
  r_ab <- cor(base,   ptime,  use = "pairwise.complete.obs")   # baseline vs pseudotime
  r_ac <- cor(method, ptime,  use = "pairwise.complete.obs")   # method   vs pseudotime
  r_bc <- cor(base,   method, use = "pairwise.complete.obs")   # overlap between predictors

  # Fisher z's
  z_ab <- 0.5 * log((1 + r_ab)/(1 - r_ab))
  z_ac <- 0.5 * log((1 + r_ac)/(1 - r_ac))

  # pooled pieces (PASS)
  rbar    <- (r_ab + r_ac)/2
  psi_bar <- r_bc * (1 - 2*rbar^2) - 0.5 * (rbar^2) * (1 - 2*rbar^2 - r_bc^2)
  sbar    <- psi_bar / (1 - rbar^2)^2

  # tiny numeric guard so SE stays real
  if (sbar >= 1) sbar <- 0.999999
  if (sbar <= -1) sbar <- -0.999999

  # pooled Steiger Z*
  se <- sqrt((2 - 2*sbar) / (n - 3))
  Z  <- (z_ac - z_ab) / se

  # one-sided p (H1 in the direction you pick)
  p_one <- if (alternative == "baseline>method") 1 - pnorm(Z) else pnorm(Z)

  data.frame(r_method = r_ac, delta_abs = abs(r_ab) - abs(r_ac),
             Z_steiger = Z, p_value = p_one, stringsAsFactors = FALSE)
}

# 2) Build per-baseline result frames in dep_list
dep_list <- vector("list", length(base))
names(dep_list) <- names(base)

for (i in seq_along(base)) {
  x_base <- base[[i]]$avg_acc
  res_i <- do.call(rbind, lapply(seq_along(method), function(j) {
    x_method <- method[[j]]$avg_acc
    out <- steiger_one_sided(x_base, x_method, n = 246, alternative = "baseline>method")
    cbind(baseline = names(base)[i], method = names(method)[j], out, stringsAsFactors = FALSE)
  }))
  dep_list[[i]] <- res_i
}

# 3) Compile to matrices: Δ|r| and one-sided p (“column > row”)
mnames <- names(method)
delta_abs_mat <- matrix(NA_real_, length(base), length(method), dimnames = list(names(base), mnames))
p_mat_one <- matrix(NA_real_, length(base), length(method), dimnames = list(names(base), mnames))

for (b in names(dep_list)) {
  df <- dep_list[[b]]
  df <- df[match(mnames, df$method), ]
  delta_abs_mat[b, ] <- df$delta_abs
  p_mat_one[b, ]     <- df$p_value
}
diag(delta_abs_mat) <- 0
diag(p_mat_one)     <- NA_real_

# 4) Heatmap with one-sided significance stars (method > baseline)
library(reshape2)
library(ggplot2)

df_delta <- melt(delta_abs_mat, varnames = c("Baseline","Method"), value.name = "delta_abs")
df_p <- melt(p_mat_one, varnames = c("Baseline","Method"), value.name = "p_one")
plot_df  <- merge(df_delta, df_p, by = c("Baseline","Method"), sort = FALSE)

plot_df$signif <- ifelse(is.na(plot_df$p_one), "",
                  ifelse(plot_df$p_one < 0.001, "***",
                  ifelse(plot_df$p_one < 0.01,  "**",
                  ifelse(plot_df$p_one < 0.05,  "*", ""))))

p <- ggplot(plot_df, aes(x = Method, y = Baseline, fill = delta_abs)) +
      geom_tile(color = "white") +
      geom_text(aes(label = signif), size = 3) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = expression(Delta*"|r|"~"(|r"[i]*"|" - "|r"[j]*"|)")) +
      coord_fixed() +
      labs(x = "Method j (column)", y = "Baseline i (row)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())

ggsave("~/Steiger_z_test.png", p, height = 5, width = 5)


## Barplot (Figure 6B)
corr_dt <- data.frame(Method = c("BART-spatial", "MERINGUE", "SpaTrack (decrease)", "SpaGene", "singleCellHaystack"), 
                      corr_eff = c(cor(BARTspatial$avg_acc, BARTspatial$spATAC_traj), cor(MERINGUE$avg_acc, MERINGUE$spATAC_traj),
                                   cor(spatrack_down$avg_acc, spatrack_down$spATAC_traj), 
                                   cor(SpaGene$avg_acc, SpaGene$spATAC_traj), cor(haystack$avg_acc, haystack$spATAC_traj)))

label_map <- c("BART-spatial" = "BART-spatial", "MERINGUE" = "MERINGUE", "SpaGene" = "SpaGene", 
               "singleCellHaystack" = "singleCell\nHaystack", "SpaTrack (decrease)" = "SpaTrack\n(decrease)")

set.seed(1)
palette <- c(sci_colors(4, effect = "contrast"), "#4C78A8", "#F58518")
names(palette) <- c("MERINGUE", "SpaGene", "SPADE", "singleCellHaystack", "BART-spatial", "SpaTrack (decrease)")

pvals <- df_p[df_p$Baseline == "BART-spatial", ]
plot_dt <- merge(corr_dt, pvals, by = "Method")

method_order <- c("BART-spatial", "MERINGUE", "SpaTrack (decrease)", "SpaGene", "singleCellHaystack")
plot_dt$Method <- factor(plot_dt$Method, levels = method_order)

plot_dt$label <- ifelse(is.na(plot_dt$p_one), "",
                      ifelse(plot_dt$p_one < 1e-3, "p < 1e-3", sprintf("p = %.3g", plot_dt$p_one)))
                      
p <- ggplot(plot_dt, aes(x = Method, y = corr_eff, fill = Method)) +
        geom_col(width = 0.7) +  
        geom_text(aes(label = label), vjust = 1.3, size = 4, family = "Arial") +
        scale_fill_manual(values = palette, guide = "none") +
        scale_x_discrete(labels = label_map) +
        labs(x = "\n", y = "Correlation coefficient") +
        theme_bw() +
        ylim(-0.55, 0) +
        geom_hline(yintercept = 0, linetype="dashed", linewidth = 0.5, color = "#404040") +
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(size = 16, face = "bold", family = "Arial", hjust = 0.5, vjust = -0.5), 
              axis.text = element_text(size = 12, family = "Arial"), 
              axis.title = element_text(size = 14, face = "bold", family = "Arial"))
ggsave("~/ptime_cor_barplot.png", p, height = 3, width = 6)


# Figure S3E
corr_dt <- data.frame(Method = c("MERINGUE", "SpaTrack (decrease)", "SpaGene", "singleCellHaystack", "BART-spatial unique"), 
                      corr_eff = c(cor(MERINGUE$avg_acc, MERINGUE$spATAC_traj),
                                   cor(spatrack_down$avg_acc, spatrack_down$spATAC_traj), 
                                   cor(SpaGene$avg_acc, SpaGene$spATAC_traj), cor(haystack$avg_acc, haystack$spATAC_traj), 
                                   cor(BARTspatial_u$avg_acc, BARTspatial_u$spATAC_traj)),
                      globalvar = c("MERINGUE", "spatrack_down", "SpaGene", "haystack", "BARTspatial_u"))

label_map <- c("BART-spatial unique" = "BART-spatial\nunique", "MERINGUE" = "MERINGUE", "SpaGene" = "SpaGene", 
               "singleCellHaystack" = "singleCell\nHaystack", "SpaTrack (decrease)" = "SpaTrack\n(decrease)")

set.seed(1)
palette <- c(sci_colors(4, effect = "contrast"), "#4C78A8", "#F58518")
names(palette) <- c("MERINGUE", "SpaGene", "SPADE", "singleCellHaystack", "BART-spatial unique", "SpaTrack (decrease)")

pvals <- df_p[df_p$Baseline == "BART-spatial unique", ]
plot_dt <- merge(corr_dt, pvals, by = "Method")

method_order <- c("BART-spatial unique", "MERINGUE", "SpaTrack (decrease)", "SpaGene", "singleCellHaystack")
plot_dt$Method <- factor(plot_dt$Method, levels = method_order)
plot_dt$label <- ifelse(is.na(plot_dt$p_one), "",
                      ifelse(plot_dt$p_one < 1e-3, "p < 1e-3", sprintf("p = %.3g", plot_dt$p_one)))
                      
p <- ggplot(plot_dt, aes(x = Method, y = corr_eff, fill = Method)) +
        geom_col(width = 0.7) +  
        geom_text(aes(label = label), vjust = 1.3, size = 4, family = "Arial") +
        scale_fill_manual(values = palette, guide = "none") +
        scale_x_discrete(labels = label_map) +
        labs(x = "\n", y = "Correlation coefficient") +
        theme_bw() +
        ylim(-0.55, 0) +
        geom_hline(yintercept = 0, linetype="dashed", linewidth = 0.5, color = "#404040") +
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(size = 16, face = "bold", family = "Arial", hjust = 0.5, vjust = -0.5), 
              axis.text = element_text(size = 12, family = "Arial"), 
              axis.title = element_text(size = 14, face = "bold", family = "Arial"))
ggsave("~/ptime_cor_barplot_bartu_base.png", p, height = 3, width = 6)



