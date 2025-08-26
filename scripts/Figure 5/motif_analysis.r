
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)

source("~/manuscript_figures.r")

spRNA_bart_result <- read.csv("~/results/spRNA/BART_results.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01) #35
spRNA_bart_result$TF <- str_to_title(tolower(spRNA_bart_result$TF))

subset <- readRDS("~/data/spRNA_subset.rds")
DefaultAssay(subset) <- "Spatial"



##### 1. Check expression of top-ranked TF coding gene
p1 <- visium_plot_gene(subset, "Sox9", 50)
p2 <- visium_plot_gene(subset, "Neurod2", 50)
p3 <- visium_plot_gene(subset, "Klf4", 50)
ggsave("Sox9_expr.png", p1, width = 3.5, height = 4, dpi = 600)
ggsave("Neurod2_expr.png", p2, width = 3.5, height = 4, dpi = 600)
ggsave("Klf4_expr.png", p3, width = 3.5, height = 4, dpi = 600)



##### 2. Check whether the coding gene belong to SVG or traj-gene
expression_matrix <- subset@assays$Spatial@counts
in_TF_name <- spRNA_bart_result[spRNA_bart_result$TF %in% rownames(expression_matrix), ] 
not_in_TF_name <- spRNA_bart_result$TF[!spRNA_bart_result$TF %in% rownames(expression_matrix)] 

traj_DEG_all <- read.csv("~/results/spRNA/traj_all_16159_genes.csv") 
traj_DEG <- traj_DEG_all %>% filter(adjusted_pvals < 0.1 & abs(correlation_rho) > 0.0921616)
traj_TF <- traj_DEG[traj_DEG$significant_features %in% spRNA_bart_result$TF, ] 

moran_DEG_all <- read.csv("~/results/spRNA/sp_all_10379_genes.csv") 
moran_DEG <- moran_DEG_all %>% filter(adjusted_pvals < 0.1 & Deviation_from_expectation > 0.07973)
sp_TF <- moran_DEG[moran_DEG$significant_features %in% spRNA_bart_result$TF, ]

overlap_TF <- intersect(traj_TF$significant_features, sp_TF$significant_features) 
traj_DEG[traj_DEG$significant_features %in% overlap_TF, ]
moran_DEG[moran_DEG$significant_features %in% overlap_TF, ]

test1 <- traj_TF[!traj_TF$significant_features %in% overlap_TF, ]
test2 <- sp_TF[!sp_TF$significant_features %in% overlap_TF, ]
moran_DEG_all[moran_DEG_all$significant_features %in% test1$significant_features, ]
traj_DEG_all[traj_DEG_all$significant_features %in% test2$significant_features, ]

in_TF <- in_TF_name[in_TF_name$TF %in% overlap_TF, ]
in_TF <- in_TF %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)

not_in_TF <- in_TF_name[!in_TF_name$TF %in% c(traj_TF$significant_features, sp_TF$significant_features), ] #14
not_in_TF <- not_in_TF %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)




##### 3. Fisher exact test to test the motif spatial distribution at UDHS-level
library(GenomicRanges)

ls <- fisher_test("KLF4", 15)
p5 <- fisher_plot(ls, "< 2.2e-16", 15)
ggsave("KLF4_motif_fisher.png", p5, width = 2.8, height = 2, dpi = 600)



##### 4. Check KLF4 motif-containing peak accessibility 
subset <- readRDS("~/spRNA_subset.rds")
DefaultAssay(subset) <- "peaks"

p6 <- motif_acc(subset, "KLF4")
ggsave("KLF4_motif_acc.png", p6, width = 4, height = 5, dpi = 600)
