
library(BARTsp)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

##### 1. Prepare input for BARTsp
E13_sp <- readRDS("~/data/spRNA_preprocess.rds")

# expression_matrix
expression_matrix <- E13_sp@assays$Spatial@counts
# meta_data
cell_metadata <- E13_sp@meta.data
cell_metadata$cell_type <- cell_metadata$predicted.id
cell_metadata$position <- ifelse(cell_metadata$cell_type == "Radial glia", "upstream", 
                            ifelse(cell_metadata$cell_type == "Postmitotic premature neurons", "downstream", NA))
# spatial_coordinates
spatial_coordinates <- data.frame(E13_sp@images$slice1@coordinates)
spatial_coordinates$x <- spatial_coordinates$imagerow
spatial_coordinates$y <- spatial_coordinates$imagecol

obj <- prepare_input(expression_matrix, cell_metadata, feature_metadata, spatial_coordinates, c("Radial glia", "Postmitotic premature neurons"))


##### 2. Pseudo-time analysis (calculate TVFs)
cds <- construct_trajectory(obj, "Radial glia")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.1250701, cor_cutoff_neg = -0.2456667) #top5%, Q1-1.5*IQR


###### 3. Spatial autocorrelation analysis (calculate SVFs)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.2373)


##### 4. Construct input for BART algorithm
genes <- get_sig_features_geneset(traj_DEG, moran_DEG)
geneset <- construct_BART_geneset_input(genes)


##### 5. Run BART algorithm
# Decreasingly-expressed genes (upstream)
bart_proj <- bart(name = "radial glia to PPN", genome = "mm10", data = geneset$upstream$feature, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

tf_interest <- c("PAX6", "SOX9", "NEUROD2", "KLF4", "FEZF2")
sig_result_TF <- results_geneset_up %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_geneset_up, c("PAX6","SOX9","NEUROD2","KLF4"), 0.05, 6)
ggsave("~/results/BART_results_upstream.pdf", p, width = 3, height = 3)

# Increasingly-expressed genes (downstream)
bart_proj <- bart(name = "radial glia to PPN", genome = "mm10", data = geneset$downstream$feature, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")

tf_interest <- c("PAX6", "SOX9", "NEUROD2", "KLF4", "FEZF2")
sig_result_TF <- results_geneset_down %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_geneset_down, c("PAX6","SOX9","NEUROD2","KLF4"), 0.05, 6)
ggsave("~/results/BART_results_downstream.pdf", p, width = 3, height = 3)



##### 6. Validate BARTsp prediction results
# Validation BART predicted TFs by pathway analysis
library(clusterProfiler)
library(org.Mm.eg.db)

results_geneset_up <- read.csv("~/results/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
results_geneset_down <- read.csv("~/results/BART_results_downstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)

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
# gene.df <- bitr(gene_set, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
geneList <- na.omit(genes$ENTREZID)
ego <- enrichGO(gene = geneList, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_result <- ego@result %>% dplyr::select(ID, Description, FoldEnrichment, p.adjust, Count) %>% filter(p.adjust < 0.05) %>% arrange(desc(Count))

# Visualization
sig_pathway <- c("growth", "central nervous system development", "developmental growth", "brain development", "negative regulation of cell differentiation", 
                 "regulation of growth", "forebrain development", "regulation of developmental growth", "cell fate commitment", "neural precursor cell proliferation", 
                 "glial cell differentiation", "cerebral cortex development", "stem cell population maintenance", "cell proliferation in forebrain", 
                 "negative regulation of neuron differentiation", "cell fate specification", "neuron fate specification", "brain development", "glial cell development",
                 "positive regulation of developmental growth", "radial glial cell differentiation", "glial cell fate specification", "negative regulation of glial cell differentiation", 
                 "cell morphogenesis involved in neuron differentiation", "neuron projection development", "positive regulation of growth", "forebrain radial glial cell differentiation", 
                 "neuron fate commitment", "forebrain cell migration")

sig_pathway_result <- ego_result %>% filter(Description %in% sig_pathway)

p <- ggplot(sig_pathway_result, aes(x = Count, y = forcats::fct_reorder(Description, Count), color = FoldEnrichment, size = -log10(p.adjust))) +
  geom_point() + 
  theme_bw() +
  labs(x = "Count", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.key.size = unit(0.3, "cm")) +
  guides(color = guide_colorbar(title = "Fold enrichment", barwidth = 1, barheight = 5),
         size = guide_legend(title = "-log10(pval_adj)"))
ggsave("~/results/downstream_TF_downstream_genes_GO_terms.pdf", p, height = 4, width = 6)


