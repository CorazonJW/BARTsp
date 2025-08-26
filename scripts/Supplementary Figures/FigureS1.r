
library(dplyr)

source("~/manucript_figures.r")


library(BARTsp)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

load("/project/zanglab_project/jw4xtu/spatial_project/tumor/data/sample_2.Robj")
object <- sample_2

##### Run BART with different genesets
## 1. Prepare input for BARTsp
expression_matrix <- object@assays$SCT@counts
cell_metadata <- object@meta.data
cell_metadata$cell_type <- as.character(cell_metadata$cluster_annotations)
spatial_coordinates <- GetTissueCoordinates(object)
colnames(spatial_coordinates) <- c("x", "y")

cell_type <- c("core", "edge", "transitory")
obj <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_type)

## 2. Spatial analysis (calculate SVFs)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}
moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.1)

## 3. Pseudo-time analysis (calculate TVFs)
cds <- construct_trajectory(obj, "core")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.15, cor_cutoff_neg = -0.1)

## 4. Construct input for BART algorithm
genes_tem <- traj_DEG
names(genes_tem) <- c("significant_features", "statistics", "cor.rho")
genes_tem <- construct_BART_geneset_input(genes_tem) 

genes_sp <- moran_DEG

genes_overlapped <- get_sig_features_geneset(traj_DEG, moran_DEG)
genes_overlapped <- construct_BART_geneset_input(genes_overlapped)

## 5. Run BART algorithm
bart_proj <- bart(name = "GRN", genome = "hg38", data = genes_tem$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")
write.csv(results_geneset_up, "BART_results_TVF.csv")

bart_proj <- bart(name = "GRN", genome = "hg38", data = genes_tem$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")
write.csv(results_geneset_up, "BART_results_SVF.csv")

bart_proj <- bart(name = "GRN", genome = "hg38", data = genes_overlapped$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")
write.csv(results_geneset_up, "BART_results_SVF.csv")




##### Plotting
oscc_up_sp_only <- read.csv("~/BART_results_SVF.csv") %>% mutate(rank = X)
oscc_up_traj_only <- read.csv("~/BART_results_TVF.csv") %>% mutate(rank = X)
oscc_up_all <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X)



# Figure S1C
oscc_up_sp_only <- scale_rank(oscc_up_sp_only, tf_interest)
oscc_up_traj_only <- scale_rank(oscc_up_traj_only, tf_interest)
oscc_up_all <- scale_rank(oscc_up_all, tf_interest)

oscc_up_sp_only$method   <- "SVFs\nonly"
oscc_up_traj_only$method <- "TVFs\nonly"
oscc_up_all$method       <- "Overlapped\nfeatures"

combined_rank <- bind_rows(
  oscc_up_sp_only     %>% select(TF, scaled_rank, TF_type, method),
  oscc_up_traj_only %>% select(TF, scaled_rank, TF_type, method), 
  oscc_up_all  %>% select(TF, scaled_rank, TF_type, method)
)
p7 <- VS_boxplot(combined_rank)
ggsave("scaled_rank.png", p7, height = 3, width = 4.5, dpi = 600)



# Figure S1D-F
object <- readRDS("/project/zanglab_project/jw4xtu/spatial_project/tumor/data/sample_2_with_ptime.Robj")

cell_df <- object@meta.data

p4 <- avg_exp_input_gene(cell_df, x, y, SVFs_expr, "SVF average\nexpression")
p5 <- avg_exp_input_gene(cell_df, x, y, TVFs_expr, "TVFs average\nexpression")
p6 <- avg_exp_input_gene(cell_df, x, y, down_gene_avg_exp, "Overlapped features\naverage expression")
ggsave("input_gene_exp_sp_only.png", p4, width = 4, height = 4.5, dpi = 600)
ggsave("input_gene_exp_traj_only.png", p5, width = 4, height = 4.5, dpi = 600)
ggsave("input_gene_exp_overlapped.png", p6, width = 4, height = 4.5, dpi = 600)



# Figure S1G-I
tf_interest <- c("GRHL3", "TCF4", "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1", "FOSL1", "STAT3", "FOS", "JUN", "JUND")

p1 <- plot_BART_results(oscc_up_sp_only, tf_interest, 0.1, 6)
p2 <- plot_BART_results(oscc_up_traj_only, tf_interest, 0.1, 6)
p3 <- plot_BART_results(oscc_up_all, tf_interest, 0.1, 6)
ggsave("BART_sp_only_prediction.png", p1, width = 3, height = 3, dpi = 600)
ggsave("BART_traj_only_prediction.png", p2, width = 3, height = 3, dpi = 600)
ggsave("BART_overlapped_prediction.png", p3, width = 3, height = 3, dpi = 600)
