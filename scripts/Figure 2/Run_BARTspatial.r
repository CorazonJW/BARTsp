
output_dir <- "~/results/"

library(BARTsp)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

source("~/manuscript_figures.r")

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

# Visualization
pallette <- c("core" = "#E5CCFF", "edge" = "#4D1A89", "transitory" = "#B266FF")

p1 <- cell_type_spatial(cell_metadata, x, y, cell_type, pallette)
ggsave("spatial_map.png", p1, width = 3.5, height = 4, dpi = 600)

p3 <- plot_gene(expression_matrix, "LCN2", cell_metadata, x, y)
p4 <- plot_gene(expression_matrix, "KRT6C", cell_metadata, x, y)
p5 <- plot_gene(expression_matrix, "COL1A2", cell_metadata, x, y)
ggsave("LCN2_expr.png", p3, width = 3.5, height = 4, dpi = 600)
ggsave("KRT6C_expr.png", p4, width = 3.5, height = 4, dpi = 600)
ggsave("COL1A2_expr.png", p5, width = 3.5, height = 4, dpi = 600)



##### 2. Pseudo-time analysis (calculate TVFs)
cds <- construct_trajectory(obj, "core")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.15, cor_cutoff_neg = -0.1)

# Visualization
pseudotime_df <- data.frame(pseudotime_values)
object <- subset(object, subset = cluster_annotations %in% cell_metadata$cell_type)
object <- AddMetaData(object, metadata = pseudotime_df$pseudotime_values, col.name = "Pseudotime")
object <- object[, is.finite(object$Pseudotime)]

cell_df <- object@cell_metadata
p2 <- ptime_spatial(cell_df, x, y, pseudotime)
ggsave("pseudotime.png", p2, width = 3.5, height = 4, dpi = 600)



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

# Visualization expression of input genes
gene_list <- rownames(obj$expression_matrix)
up_sig_genes <- geneset$up_gene$significant_features
down_sig_genes <- geneset$down_gene$significant_features

up_expr <- Matrix::colMeans(obj$expression_matrix[up_sig_genes, , drop = FALSE])
down_expr <- Matrix::colMeans(obj$expression_matrix[down_sig_genes, , drop = FALSE])

cell_df$up_gene_avg_exp <- up_expr
cell_df$down_gene_avg_exp <- down_expr

p6 <- avg_exp_input_gene(cell_df, x, y, up_gene_avg_exp)
p7 <- avg_exp_input_gene(cell_df, x, y, down_gene_avg_exp)
ggsave("input_gene_exp_up.png", p6, width = 3.5, height = 4, dpi = 600)
ggsave("input_gene_exp_down.png", p7, width = 3.5, height = 4, dpi = 600)



##### 5. Run BART algorithm
# Decreasingly-expressed genes (upstream - core)
bart_proj <- bart(name = "oscc", genome = "hg38", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")
write.csv(results_geneset_up, "~/BART_results_upstream.csv")

tf_interest <- c("GRHL3", "TCF4", 
                 "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                 "FOS", "JUN", "JUND", "TP53")

print(results_geneset_up[results_geneset_up$TF %in% tf_interest, ])
print(sig_result_up[sig_result_up$TF %in% tf_interest, ])

p <- plot_BART_results(results_geneset_up, tf_interest, 0.1, 6)
ggsave("BART_results_upstream.png", p, width = 3, height = 3, dpi = 600)

# Increasingly-expressed genes (downstream)
bart_proj <- bart(name = "oscc", genome = "hg38", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")
write.csv(results_geneset_down, "~/BART_results_downstream.csv")

tf_interest <- c("GRHL3", "TCF4", 
                 "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                 "FOS", "JUN", "JUND", "TP53")

print(results_geneset_down[results_geneset_down$TF %in% tf_interest, ])
print(sig_result_down[sig_result_down$TF %in% tf_interest, ])

p <- plot_BART_results(results_geneset_down, tf_interest, 0.1, 6)
ggsave("BART_results_downstream.png", p, width = 3, height = 3, dpi = 600)



##### 6. Validate BARTsp prediction results by pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

results_geneset_up <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
# results_geneset_down <- read.csv("~/BART_results_downstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
results_geneset_up$TF <- str_to_title(tolower(results_geneset_up$TF))

mouse_TF_target <- read.table("~/TRRUST/trrust_rawdata.human.tsv")
colnames(mouse_TF_target) <- c("TF", "target", "relationship", "ref")
mouse_TF_target <- mouse_TF_target[which(mouse_TF_target$relationship %in% c("Activation")), -4]
TRRUST_TF <- mouse_TF_target %>% filter(TF %in% results_geneset_up$TF)
nrow(TRRUST_TF)

gene_set <- TRRUST_TF$target
genes <- select(org.Hs.eg.db, keys = gene_set, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
geneList <- na.omit(genes$ENTREZID)
ego <- enrichGO(gene = geneList, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_result <- ego@result %>% dplyr::select(ID, Description, FoldEnrichment, p.adjust, Count) %>% filter(p.adjust < 0.05) %>% arrange(desc(Count))

up_TF_ego_result <- data.frame(ego_result)

# Visualization
library(dplyr)
library(ggplot2)
library(forcats)

ego_sig_result_1 <- up_TF_ego_result %>% filter(grepl("STAT|Wnt|MAPK|wound healing", Description, ignore.case = FALSE))
ego_sig_result_2 <- down_TF_ego_result %>% filter(grepl("STAT|Wnt|MAPK|wound healing", Description, ignore.case = FALSE))
sig_pathway_result_1 <- ego_sig_result_1 %>% filter(Description %in% ego_sig_result_1$Description)
sig_pathway_result_2 <- ego_sig_result_2 %>% filter(Description %in% ego_sig_result_2$Description)

sig_pathway_result_1$Group <- "Upstream"
sig_pathway_result_2$Group <- "Downstream"
sig_pathway_result_1$Group <- factor(sig_pathway_result_1$Group, levels = c("Upstream", "Downstream"))
sig_pathway_result_2$Group <- factor(sig_pathway_result_2$Group, levels = c("Upstream", "Downstream"))

combined_result <- bind_rows(head(sig_pathway_result_1, 20), head(sig_pathway_result_2, 20))

p <- ggplot(combined_result, aes(x = Count, y = forcats::fct_reorder(Description, Count), 
                                 color = FoldEnrichment, size = -log10(p.adjust))) +
  geom_point() + 
  theme_bw() +
  labs(x = "Count", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  theme(legend.position = "right", 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 9), 
        legend.key.size = unit(0.3, "cm"), 
        strip.text = element_text(size = 10, face = "bold")) +
  guides(color = guide_colorbar(title = "Fold enrichment", barwidth = 1, barheight = 5),
         size = guide_legend(title = "-log10(pval_adj)"))
ggsave("~/TF_downstream_genes_GO_terms.pdf", p, height = 4, width = 6)

