
output_dir <- "~/results/"

library(BARTsp)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

##### 1. Prepare input for BARTsp
load("~/data/sample_2.Robj")
object <- sample_2

# expression_matrix
expression_matrix <- object@assays$SCT@counts
# meta_data
cell_metadata <- object@meta.data
cell_metadata$cell_type <- as.character(cell_metadata$cluster_annotations)
# spatial_coordinates
spatial_coordinates <- GetTissueCoordinates(object)
colnames(spatial_coordinates) <- c("x", "y")

cell_type <- c("core", "edge", "transitory")
obj <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_type)

##### 2. Pseudo-time analysis (calculate TVFs)
cds <- construct_trajectory(obj, "core")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.15, cor_cutoff_neg = -0.1)

# Visualization
pseudotime_df <- data.frame(pseudotime_values)
object <- subset(object, subset = cluster_annotations %in% cell_metadata$cell_type)
object <- AddMetaData(object, metadata = pseudotime_df$pseudotime_values, col.name = "Pseudotime")
object <- object[, is.finite(object$Pseudotime)]

cell_df <- obj$cell_metadata
cell_df$x <- obj$spatial_coordinates[, 2]
cell_df$y <- obj$spatial_coordinates[, 3]
cell_df$pseudotime <- pseudotime_df$pseudotime_values

pallette <- c("core" = "#E41A1C", "edge" = "#4DAF4A", "transitory" = "#FFFF33")

p1 <- ggplot(cell_df, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = pallette) +
  labs(color = "Region", x = "X Coordinate", y = "Y Coordinate") +
  theme_bw() +  
  theme(legend.position = "top", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.text = element_blank(), axis.ticks = element_blank())
p2 <- ggplot(cell_df, aes(x = x, y = y, color = pseudotime_values)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  labs(color = "Pseudotime", x = "X Coordinate", y = "Y Coordinate") +
  theme_bw() +  
  theme(legend.position = "top", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.text = element_blank(), axis.ticks = element_blank())
p <- p1+p2
ggsave("~/results/pseudotime.pdf", p, width = 6, height = 4)


###### 3. Spatial autocorrelation analysis (calculate SVFs)
obj$expression_matrix <- as.matrix(obj$expression_matrix)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.1)


##### 4. Construct input for BART algorithm
genes <- get_sig_features_geneset(traj_DEG, moran_DEG)
geneset <- construct_BART_geneset_input(genes)

# Check expression of input genes
gene_list <- rownames(obj$expression_matrix)
up_sig_genes <- geneset$up_gene$significant_features
down_sig_genes <- geneset$down_gene$significant_features

up_expr <- Matrix::colMeans(obj$expression_matrix[up_sig_genes, , drop = FALSE])
down_expr <- Matrix::colMeans(obj$expression_matrix[down_sig_genes, , drop = FALSE])

cell_df$up_gene_avg_exp <- up_expr
cell_df$down_gene_avg_exp <- down_expr

# Up-regulated genes
p1 <- ggplot(cell_df, aes(x = x, y = y, color = up_gene_avg_exp)) +
        geom_point(size = 1.5) +
        scale_color_viridis_c() +
        theme_bw() +
        theme(legend.position = "top")
# Down-regulated genes
p2 <- ggplot(cell_df, aes(x = x, y = y, color = down_gene_avg_exp)) +
        geom_point(size = 1.5) +
        scale_color_viridis_c() +
        theme_bw() +
        theme(legend.position = "top")
p <- p1 + p2
ggsave("~/results/input_gene_exp.pdf", p, width = 7.5, height = 4.5)



##### 5. Run BART algorithm
# Decreasingly-expressed genes (upstream - core)
bart_proj <- bart(name = "oscc", genome = "hg38", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

tf_interest <- c("GRHL3", "TCF4", 
                 "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                 "FOS", "JUN", "JUND", "TP53")

print(results_geneset_up[results_geneset_up$TF %in% tf_interest, ])
print(sig_result_up[sig_result_up$TF %in% tf_interest, ])

p <- plot_BART_results(results_geneset_up, tf_interest, 0.1, 6)
ggsave(paste0(output_dir, "BART_results_upstream.pdf"), p, width = 3, height = 3)

# Increasingly-expressed genes (downstream)
bart_proj <- bart(name = "oscc", genome = "hg38", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")


tf_interest <- c("GRHL3", "TCF4", 
                 "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                 "FOS", "JUN", "JUND", "TP53")

print(results_geneset_down[results_geneset_down$TF %in% tf_interest, ])
print(sig_result_down[sig_result_down$TF %in% tf_interest, ])

p <- plot_BART_results(results_geneset_down, tf_interest, 0.1, 6)
ggsave(paste0(output_dir, "BART_results_downstream.pdf"), p, width = 3, height = 3)





##### 6. Validate BARTsp prediction results
# Check difference in ranks
library(stringr)
library(ggsignif)

results_geneset_up <- read.csv("~/results/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
results_geneset_up$TF_type <- ifelse(results_geneset_up$TF %in% tf_interest, "functional", "not_functional")
results_geneset_up$TF_type <- factor(results_geneset_up$TF_type, levels = c("functional", "not_functional"))
results_geneset_up$rank <- as.numeric(results_geneset_up$rank)
results_geneset_up$scaled_rank <- round((max(results_geneset_up$rank) - results_geneset_up$rank) / (max(results_geneset_up$rank) - min(results_geneset_up$rank)), 3)

wilcox_result <- wilcox.test(results_geneset_up$scaled_rank ~ results_geneset_up$TF_type, data = results_geneset_up, exact = FALSE)
p_val <- round(wilcox_result$p.value, 2)

if (p_val <= 0.1) {
  annotation_label <- if (p_val < 0.01) "***"
                    else if (p_val < 0.05) "**"
                    else "*"
} else {
  annotation_label <- "NS"
}

p <- ggplot(results_geneset_up, aes(x = TF_type, y = scaled_rank, fill = TF_type)) +
      geom_boxplot() +
      geom_signif(comparisons = list(c("functional", "not_functional")), annotations = annotation_label, map_signif_level = TRUE, textsize = 5, vjust = 0.5) +
      theme_bw() +
      labs(title = "", x = "TF type", y = "Scaled rank") +
      theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 12))
ggsave("~/results/BART_results_upstream_TF_rank.pdf", p, width = 3.5, height = 3.5)


# Check expression of top 10 TF from up and down sets
library(gridExtra)

common_cells <- intersect(colnames(expression_matrix), rownames(cell_df))
expression_matrix <- expression_matrix[, common_cells]
cell_df <- cell_df[common_cells, ]

plot_gene <- function(gene_name) {
  expr <- expression_matrix[gene_name, ]
  upper_lim <- quantile(expr, 0.99, na.rm = TRUE)
  df <- cell_df
  df$expr <- expr
  
  p <- ggplot(df, aes(x = x, y = y, color = expr)) +
    geom_point(size = 1.0) +
    scale_color_gradientn(colors = c("#FFF6F6", "#FFCCCC", "#FF9999", "#FF3333", "#FF0000"), limits = c(0, upper_lim), 
                          oob = scales::squish, name = gene_name) +
    theme_bw() +
    theme(legend.position = "top", legend.key.height = unit(0.3, "cm"),legend.key.width = unit(0.7, "cm"),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 8),
          plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  return(p)
}

top_up_genes <- sig_result_up$TF[1:10]
top_down_genes <- sig_result_down$TF[1:10]
top_genes <- c(top_up_genes, top_down_genes)
top_genes <- intersect(top_genes, rownames(expression_matrix))

# tf_interest <- c("GRHL3", "TCF4", "TP63", "GRHL2", "SOX2", "KLF4", "TP73", "SNAI2", "ZEB1",  "FOSL1", "STAT3", "FOS", "JUN", "JUND")
# tf_interest <- c("COL1A1", "FN1", "COL1A2", "TIMP1", "COL6A2", 
#                  "KRT6C", "LYPD3", "KRTDAP", "KRT6B", "SLPI", 
#                  "SPRR2D", "SPRR2E", "DEFB4A", "SPRR2A", "LCN2", 
#                  "CLDN4", "CSTA", "SPRR1B", "LAMC2", "ITGA5")

# top_genes <- intersect(tf_interest, rownames(expression_matrix))

pdf("~/results/predicted_TF_exp_2.pdf", width = 4.5, height = 5.5)
n <- length(top_genes)
for (i in seq(1, n, by = 4)) {
  plots <- lapply(top_genes[i:min(i+3, n)], plot_gene)
  grid.arrange(grobs = plots, ncol = 2, nrow = 2)
}
dev.off()