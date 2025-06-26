
##### Validation BARTsp and SpaTrack predicted TFs by pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)

BARTsp_results <- read.csv("~/results/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.1)
spaTrack_results <- read.csv("~/results/spatrack_TF_prediction_spatrack_TF_decrease.csv") %>% mutate(rank = X) 

dataset <- spaTrack_results

mouse_TF_target <- read.table("~/data/TRRUST/trrust_rawdata.human.tsv")
colnames(mouse_TF_target) <- c("TF", "target", "relationship", "ref")
mouse_TF_target <- mouse_TF_target[which(mouse_TF_target$relationship %in% c("Activation")), -4]
TRRUST_TF <- mouse_TF_target %>% filter(TF %in% dataset$TF)
nrow(TRRUST_TF)

gene_set <- TRRUST_TF$target
genes <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_set, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
geneList <- na.omit(genes$ENTREZID)
ego <- enrichGO(gene = geneList, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_result <- ego@result %>% dplyr::select(ID, Description, FoldEnrichment, p.adjust, Count) %>% filter(p.adjust < 0.05) %>% arrange(desc(Count))

# Visualization
BARTsp_ego_results <- read.table("~/results/BARTsp_TF_downstream_genes_GO_terms.txt") %>% filter(Count > 10)
SpaTrack_ego_results <- read.table("~/results/SpaTrack_TF_downstream_genes_GO_terms.txt") # %>% filter(Count > 10)

ego_sig_result_1 <- SpaTrack_ego_results %>% filter(grepl("STAT|Wnt|MAPK|wound healing", Description, ignore.case = FALSE))
ego_sig_result_2 <- BARTsp_ego_results %>% filter(grepl("STAT|Wnt|MAPK|wound healing", Description, ignore.case = FALSE))

sig_pathway_result_1 <- ego_sig_result_1 %>% filter(Description %in% ego_sig_result_1$Description)
sig_pathway_result_2 <- ego_sig_result_2 %>% filter(Description %in% ego_sig_result_2$Description)

library(dplyr)
library(ggplot2)
library(forcats)

# Label and combine data
sig_pathway_result_1$Group <- "SpaTrack"
sig_pathway_result_2$Group <- "BARTsp"

# Set factor levels to enforce order: Up first, then Down
sig_pathway_result_1$Group <- factor(sig_pathway_result_1$Group, levels = c("SpaTrack", "BARTsp"))
sig_pathway_result_2$Group <- factor(sig_pathway_result_2$Group, levels = c("SpaTrack", "BARTsp"))

# Combine and select top 20 from each
combined_result <- bind_rows(head(sig_pathway_result_1, 20), head(sig_pathway_result_2, 20))
combined_result$FC <- pmin(combined_result$FoldEnrichment, 10)

# Plot with facet_grid for vertical layout
p <- ggplot(combined_result, aes(x = Count, y = forcats::fct_reorder(Description, Count), color = FC, size = -log10(p.adjust))) +
  geom_point() + 
  theme_bw() +
  labs(x = "Count", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 10), legend.key.size = unit(0.2, "cm"), strip.text = element_text(size = 10, face = "bold")) +
  guides(color = guide_colorbar(title = "Fold enrichment", barwidth = 5, barheight = 1), size = guide_legend(title = "-log10(pval_adj)"))
ggsave("~/results/TF_downstream_genes_GO_terms.pdf", p, height = 5.5, width = 4.5)
