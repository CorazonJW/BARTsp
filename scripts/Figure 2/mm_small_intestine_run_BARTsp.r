

library(Seurat)
library(future)
plan("multisession", workers = 10)

output_dir <- "~/results"

library(BARTsp)

library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

##### 1. Prepare input for BARTsp
object <- readRDS("~/subset_object_with_ct_label.RDS")

# expression_matrix
expression_matrix <- object@assays$Spatial.008um@layers$counts
colnames(expression_matrix) <- rownames(object@assays$Spatial.008um@cells)
rownames(expression_matrix) <- rownames(object@assays$Spatial.008um@features)
# meta_data
cell_metadata <- object@meta.data
cell_metadata$Cell_type <- cell_metadata$Cell_type
cell_metadata$cell_type <- ifelse(cell_metadata$Cell_type %in% c("Enterocyte.Progenitor", "Enterocyte.Progenitor.Late", "Enterocyte.Progenitor.Early"), "Enterocyte_Progenitor",
                                        ifelse(cell_metadata$Cell_type == "Enterocyte.Immature.Distal", "Enterocyte_Immature", 
                                            ifelse(cell_metadata$Cell_type == "Enterocyte.Mature.Distal", "Enterocyte_Mature", NA)))
# spatial_coordinates
spatial_coordinates <- GetTissueCoordinates(object)

obj <- prepare_input(expression_matrix, cell_metadata, feature_metadata, spatial_coordinates, c("Enterocyte_Progenitor", "Enterocyte_Immature", "Enterocyte_Mature"))



##### 2. Pseudo-time analysis (calculate TVFs)
cds <- construct_trajectory(obj, "Enterocyte_Progenitor")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.05, cor_cutoff_neg = -0.05)

# Visualization
pseudotime_df <- data.frame(pseudotime_values)
object <- subset(object, subset = enterocyte_type %in% cell_metadata$cell_type)
object <- AddMetaData(object, metadata = pseudotime_df$pseudotime_values, col.name = "Pseudotime")
object <- object[, is.finite(object$Pseudotime)]

pallette <- c("Enterocyte_Immature" = "#E41A1C", "Enterocyte_Mature" = "#4DAF4A", "Enterocyte_Progenitor" = "#FFFF33")

Idents(object) <- object$enterocyte_type
p1 <- SpatialFeaturePlot(object, "Pseudotime", pt.size.factor = 18, image.alpha = 0.7, alpha = c(0.8, 0.8)) +
        theme(legend.key.size = unit(0.8, "cm"))
p2 <- SpatialDimPlot(object, label = F, pt.size.factor = 18, image.alpha = 0.7, alpha = c(0.8, 0.8), cols = pallette) + 
        labs(fill = "Cell Type") + theme(legend.text = element_text(size = 12))
p <- p1 + p2
ggsave(paste0(output_dir, "enterocyte_pseudotime.pdf"), p, width = 8.5, height = 4)



###### 3. Spatial autocorrelation analysis (calculate SVFs)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.01)




##### 4. Construct input for BART algorithm
genes <- get_sig_features_geneset(traj_DEG, moran_DEG)
geneset <- construct_BART_geneset_input(genes) 

# Check expression of input genes
g <- head(genes$significant_features, 50)
for (i in g) {
     if (!i %in% rownames(object)) {
        message(paste("Skipping gene:", i, " - not found"))
        next
     }

     expr_values <- FetchData(object, vars = i, slot = "data")[[1]]
     upper_lim <- quantile(expr_values, 0.99, na.rm = TRUE)

     p <- SpatialFeaturePlot(object, features = i, pt.size.factor = 18, slot = "data") + 
            scale_fill_gradientn(colors = c("#FFF6F6", "#FFCCCC", "#FF9999", "#FF3333", "#FF0000"), limits = c(0, upper_lim), oob = scales::squish) + 
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), 
                  legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.5, "cm"))
     ggsave(paste0(output_dir, "/input_gene_exp/", i, "_exp.pdf"), p, width = 5, height = 4)
     
     message("done: ", i)
}


# Check pathways enriched of input genes to BARTsp
library(clusterProfiler)
library(org.Mm.eg.db)

gene_set <- genes$significant_features
genes <- select(org.Mm.eg.db, keys = gene_set, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
geneList <- na.omit(genes$ENTREZID)
ego <- enrichGO(gene = geneList, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
ego_result <- ego@result %>% dplyr::select(ID, Description, FoldEnrichment, p.adjust, Count) %>% filter(p.adjust < 0.1) %>% arrange(desc(Count))

# Visualization
ego_sig_result <- ego_result %>% filter(grepl("intestine|intestinal|epithelial|enterocyte", Description, ignore.case = TRUE))
sig_pathway_result <- ego_result %>% filter(Description %in% ego_sig_result$Description)

p <- ggplot(head(sig_pathway_result, 20), aes(x = Count, y = forcats::fct_reorder(Description, Count), color = FoldEnrichment, size = -log10(p.adjust))) +
  geom_point() + 
  theme_bw() +
  labs(x = "Count", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.key.size = unit(0.3, "cm")) +
  guides(color = guide_colorbar(title = "Fold enrichment", barwidth = 1, barheight = 5),
         size = guide_legend(title = "-log10(pval_adj)"))
ggsave(paste0(output_dir, "/input_genes_GO_terms.pdf"), p, height = 4, width = 6)





##### 5. Run BART algorithm
# Decreasingly-expressed genes (upstream)
bart_proj <- bart(name = "enterocyte", genome = "mm10", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

results_geneset_up <- results_geneset_up %>% mutate(rank = rownames(.))
sig_result_up <- results_geneset_up %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.05)
nrow(sig_result_up)

tf_interest <- c("HNF4G", "HNF4A", "GATA6", "HNF1B", "MAF", "CDX2", "MAFB", "GATA4", "HES1")
sig_result_TF <- results_geneset_up %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_geneset_up, tf_interest, 0.05, 6)
ggsave(paste0(output_dir, "~/BART_results_upstream.pdf"), p, width = 3, height = 3)

# Increasingly-expressed genes (downstream)
bart_proj <- bart(name = "enterocyte", genome = "mm10", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")

results_geneset_down <- results_geneset_down %>% mutate(rank = rownames(.))
sig_result_down <- results_geneset_down %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.05)
nrow(sig_result_down) 

tf_interest <- c("HNF4G", "HNF4A", "GATA6", "HNF1B", "MAF", "CDX2", "MAFB", "GATA4", "HES1")
sig_result_TF <- results_geneset_down %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_geneset_down, tf_interest, 0.05, 6)
ggsave(paste0(output_dir, "~/BART_results_downstream.pdf"), p, width = 3, height = 3)





##### 6. Validate BARTsp prediction results
sig_result_down <- read.csv(paste0(output_dir, "~/BART_results_downstream.csv")) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01)
sig_result_up <- read.csv(paste0(output_dir, "~/BART_results_upstream.csv")) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01)

genes_of_interest_down <- str_to_title(tolower(sig_result_down$TF)) 
genes_of_interest_up <- str_to_title(tolower(sig_result_up$TF))

entero_object <- subset(object, subset = Cell_type %in% c("Enterocyte.Progenitor", "Enterocyte.Progenitor.Late", "Enterocyte.Progenitor.Early", 
                                                          "Enterocyte.Immature.Distal", "Enterocyte.Mature.Distal"))

# Check expression of predicted TF
g <- c(genes_of_interest_down[genes_of_interest_down %in% genes$significant_features], genes_of_interest_up[genes_of_interest_up %in% genes$significant_features])

for (i in g) {
     if (!i %in% rownames(entero_object)) {
        message(paste("Skipping gene:", i, " - not found"))
        next
     }

     expr_values <- FetchData(entero_object, vars = i, slot = "data")[[1]]
     upper_lim <- quantile(expr_values, 0.999, na.rm = TRUE)

     p <- SpatialFeaturePlot(entero_object, features = i, pt.size.factor = 18, slot = "data") + 
            scale_fill_gradientn(colors = c("#FFF6F6", "#FFCCCC", "#FF9999", "#FF3333", "#FF0000"), limits = c(0, upper_lim), oob = scales::squish) + 
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), 
                  legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.5, "cm"))
     ggsave(paste0(output_dir, "/TF_conding_gene_exp/", i, "_exp.pdf"), p, width = 5, height = 4)
     
     message("done: ", i)
}



# Check pathways enriched of downstream targets of BARTsp-predicted targets
library(clusterProfiler)
library(org.Mm.eg.db)

results_geneset_up <- read.csv(paste0(output_dir, "/BART_results_upstream.csv")) %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
results_geneset_down <- read.csv(paste0(output_dir, "/BART_results_downstream.csv")) %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)

spRNA_bart_result <- results_geneset_up
spRNA_bart_result <- results_geneset_down

library(stringr)
spRNA_bart_result$TF <- str_to_title(tolower(spRNA_bart_result$TF))

mouse_TF_target <- read.table("~/data/trrust_rawdata.mouse.tsv")
colnames(mouse_TF_target) <- c("TF", "target", "relationship", "ref")
mouse_TF_target <- mouse_TF_target[which(mouse_TF_target$relationship %in% c("Activation")), -4]
TRRUST_TF <- mouse_TF_target %>% filter(TF %in% spRNA_bart_result$TF)
nrow(TRRUST_TF)

gene_set <- TRRUST_TF$target
genes <- select(org.Mm.eg.db, keys = gene_set, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
geneList <- na.omit(genes$ENTREZID)
ego <- enrichGO(gene = geneList, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_result <- ego@result %>% dplyr::select(ID, Description, FoldEnrichment, p.adjust, Count) %>% filter(p.adjust < 0.05) %>% arrange(desc(Count))

dt <- data.frame(ego_result)
write.table(dt, paste0(output_dir, "/TF_upstream_genes_GO_terms.txt"))
write.table(dt, paste0(output_dir, "/TF_downstream_genes_GO_terms.txt"))

# Visualization
up_TF_ego_result <- read.table(paste0(output_dir, "/TF_upstream_genes_GO_terms.txt"))
down_TF_ego_result <- read.table(paste0(output_dir, "/TF_downstream_genes_GO_terms.txt"))

ego_sig_result_1 <- up_TF_ego_result %>% filter(grepl("intestine|intestinal|enterocyte", Description, ignore.case = TRUE))
ego_sig_result_2 <- down_TF_ego_result %>% filter(grepl("intestine|intestinal|enterocyte", Description, ignore.case = TRUE))

sig_pathway_result_1 <- ego_sig_result_1 %>% filter(Description %in% ego_sig_result_1$Description)
sig_pathway_result_2 <- ego_sig_result_2 %>% filter(Description %in% ego_sig_result_2$Description)

library(dplyr)
library(ggplot2)
library(forcats)

sig_pathway_result_1$Group <- "Upstream"
sig_pathway_result_2$Group <- "Downstream"
sig_pathway_result_1$Group <- factor(sig_pathway_result_1$Group, levels = c("Upstream", "Downstream"))
sig_pathway_result_2$Group <- factor(sig_pathway_result_2$Group, levels = c("Upstream", "Downstream"))

combined_result <- bind_rows(head(sig_pathway_result_1, 20), head(sig_pathway_result_2, 20))

p <- ggplot(combined_result, aes(x = Count, y = forcats::fct_reorder(Description, Count), color = FoldEnrichment, size = -log10(p.adjust))) +
  geom_point() + 
  theme_bw() +
  labs(x = "Count", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  theme(legend.position = "right", axis.text.y = element_text(size = 10), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 9), legend.key.size = unit(0.3, "cm"), strip.text = element_text(size = 10, face = "bold")) +
  guides(color = guide_colorbar(title = "Fold enrichment", barwidth = 1, barheight = 5), size = guide_legend(title = "-log10(pval_adj)"))
ggsave(paste0(output_dir, "/TF_downstream_genes_GO_terms.pdf"), p, height = 4, width = 6)

