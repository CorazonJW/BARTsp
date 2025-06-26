
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
expression_matrix <- E13_sp@assays$peaks@counts
# cell_metadata
cell_metadata <- E13_sp@meta.data
cell_metadata$cell_type <- cell_metadata$predicted.id
# spatial_coordinates
spatial_coordinates <- data.frame(E13_sp@images$slice1@coordinates)
spatial_coordinates$x <- spatial_coordinates$imagerow
spatial_coordinates$y <- spatial_coordinates$imagecol

obj <- prepare_input(expression_matrix, cell_metadata, feature_metadata, spatial_coordinates, c("Radial glia", "Postmitotic premature neurons"))


##### 2. Pseudo-time analysis (calculate TVFs)
pseudotime_values <- E13_sp$spATAC_traj
traj_DAR <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, cor_cutoff_pos = 0.1993, cor_cutoff_neg = -0.2216) #457, minimum, mean

###### 3. Spatial autocorrelation analysis (calculate SVFs)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DAR <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.099)


##### 4. Construct input for BART algorithm
region <- get_sig_features_region(obj, traj_DAR, moran_DAR)
regions <- construct_BART_region_input(region, "trajectory_up", "trajectory_down")

# Check expression of input peaks
peak_list <- rownames(obj$expression_matrix)
down_region <- region$down_regions$significant_features
down_expr <- Matrix::colMeans(obj$expression_matrix[down_region, , drop = FALSE])

cell_df <- obj$cell_metadata
cell_df$x <- obj$spatial_coordinates[, 2]
cell_df$y <- obj$spatial_coordinates[, 3]
cell_df$down_peak_avg_acc <- down_expr

# Down-regulated genes
object <- AddMetaData(subset, metadata = cell_df$down_peak_avg_acc, col.name = "Average accessibility (Dec)")

p <- SpatialFeaturePlot(object, features = "Average accessibility (Dec)", pt.size.factor = 50, slot = "data") + 
            scale_fill_gradientn(colors = c("#FFF6F6", "#FFF6F6", "#FFCCCC", "#FF9999", "#FF3333", "#FF0000"), , limits = c(0, 0.8), oob = scales::squish) + 
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), 
                  legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.5, "cm"))
ggsave("~/results/input_peak_exp.pdf", p, width = 4, height = 4)




##### 5. Run BART algorithm
# Increasingly-expressed regions (downstream)
bart_proj <- bart(name = "radial glia to PPN", genome = "mm10", data = regions$up_region, type = "region")
bart_proj <- run_BART(bart_proj, type = "region")
results_region_up <- get_BART_results(bart_proj, "region")

tf_interest <- c("PAX6", "SOX9", "NEUROD2", "KLF4")
sig_result_TF <- results_region_up %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_region_up, c("PAX6","SOX9","NEUROD2","KLF4"), 0.05, 6)
ggsave("~/results/BART_results_downstream.pdf", p, width = 3, height = 3)


# Decreasingly-expressed regions (upstream)
bart_proj <- bart(name = "radial glia to PPN", genome = "mm10", data = regions$down_region, type = "region")
bart_proj <- run_BART(bart_proj, type = "region")
results_region_down <- get_BART_results(bart_proj, "region")

tf_interest <- c("PAX6", "SOX9", "NEUROD2", "KLF4")
sig_result_TF <- results_region_down %>% filter(TF %in% tf_interest)
print(sig_result_TF)

p <- plot_BART_results(results_region_down, c("PAX6","SOX9","NEUROD2","KLF4"), 0.05, 6)
ggsave("~/results/BART_results_upstream.pdf", p, width = 3, height = 3)





##### 6. Validate BARTsp prediction results
# Validation input peaks by pathway analysis
library(GenomicRanges)
library(rGREAT)
input_peaks <- read.csv("~/results/overlap_regionset_BART_input.csv")
input_peaks <- input_peaks$X

# convert to gr object
gr_significant <- GRanges(
  seqnames = gsub("-.*", "", input_peaks), 
  ranges = IRanges(start = as.numeric(gsub(".*-(\\d+)-\\d+", "\\1", input_peaks)),
                   end = as.numeric(gsub(".*-\\d+-(\\d+)", "\\1", input_peaks)))
)

# run GREAT
res_BP <- great(gr_significant, "GO:BP", "TxDb.Mmusculus.UCSC.mm10.knownGene")
tb_BP <- getEnrichmentTable(res_BP)

BP_result <- tb_BP %>% dplyr::select(description, genome_fraction, observed_region_hits, fold_enrichment, p_value, p_adjust) %>% filter(fold_enrichment >= 1 & p_adjust < 0.05)

sig_pathway <- c("central nervous system development", "brain development", "negative regulation of cell differentiation", 
                 "forebrain development", "cell fate commitment", "neural precursor cell proliferation", 
                 "glial cell differentiation","cell proliferation in forebrain", 
                 "negative regulation of neuron differentiation", "cell fate specification", "neuron fate specification", "glial cell development",
                 "radial glial cell differentiation", "glial cell fate specification", "negative regulation of glial cell differentiation", 
                 "cell morphogenesis involved in neuron differentiation", "neuron projection development", "forebrain radial glial cell differentiation", 
                 "neuron fate commitment", "forebrain cell migration", "glial cell fate commitment")

sig_pathway_result <- BP_result %>% filter(description %in% sig_pathway)
nrow(sig_pathway_result)

# Visualization
library(ggplot2)
library(forcats)
p <- ggplot(sig_pathway_result, aes(x = observed_region_hits, y = fct_reorder(description, observed_region_hits), color = fold_enrichment, size = observed_region_hits)) +
  geom_point() + 
  theme_bw() +
  labs(x = "Observed Region Hits", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", axis.text.y = element_text(size = 11)) +
  guides(color = guide_colorbar(title = "Fold enrichment"), size = guide_legend(title = "-log10(pval_adj)"))
ggsave("~/results/BART_input_regions_GO_result.pdf", p, height = 4, width = 6.5)


# Validation BART predicted TFs by pathway analysis
library(clusterProfiler)
library(org.Mm.eg.db)

spATAC_bart_result <- read.csv("~/results/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1) #119

library(stringr)
spATAC_bart_result$TF <- str_to_title(tolower(spATAC_bart_result$TF))

mouse_TF_target <- read.table("~/data/trrust_rawdata.mouse.tsv")
colnames(mouse_TF_target) <- c("TF", "target", "relationship", "ref")
mouse_TF_target <- mouse_TF_target[which(mouse_TF_target$relationship %in% c("Activation")), -4]
TRRUST_TF <- mouse_TF_target %>% filter(TF %in% spATAC_bart_result$TF)
nrow(TRRUST_TF)

gene_set <- TRRUST_TF$target
genes <- select(org.Mm.eg.db, keys = gene_set, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
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
ggsave("~/results/upstream_TF_downstream_genes_GO_terms.pdf", p, height = 4, width = 6)
