
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)

##### 1. Load prediction results from spRNA and spATAC data
spRNA_bart_result_up <- read.csv("~/results/spRNA/BART_results_upstream.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
spATAC_bart_result_up <- read.csv("~/results/spATAC/BART_results_upstream.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)

spRNA_sig <- spRNA_bart_result_up %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.05)
spATAC_sig <- spATAC_bart_result_up %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.05)

common <- intersect(spRNA_sig$TF, spATAC_sig$TF)

##### 2. Scale TF rank into range 0-1
spRNA_bart_result_up$scaled_rank <- round((max(spRNA_bart_result_up$rank) - spRNA_bart_result_up$rank) / (max(spRNA_bart_result_up$rank) - min(spRNA_bart_result_up$rank)), 3)
spATAC_bart_result_up$scaled_rank <- round((max(spATAC_bart_result_up$rank) - spATAC_bart_result_up$rank) / (max(spATAC_bart_result_up$rank) - min(spATAC_bart_result_up$rank)), 3)

##### 3. Check correlation between scaled ranks
sig_dt <- merge(spRNA_bart_result_up, spATAC_bart_result_up, by = "TF") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue.x, scaled_rank.x, rank_avg_z_p_a_irwinhall_pvalue.y, scaled_rank.y)
colnames(sig_dt) <- c("TF", "gene_pval", "gene_scaled_rank", "peak_pval", "peak_scaled_rank")

dt <- sig_dt
colnames(dt) <- c("TF", "pval1", "scale_rank_1", "pval2", "scale_rank_2")
dt <- dt %>% mutate(overlap_TF = ifelse(pval1 < 0.05 & pval2 < 0.05, "True", "False"))

correlation_res <- cor.test(dt$scale_rank_1, dt$scale_rank_2, method = "spearman")
correlation <- round(correlation_res$estimate, 3)
print(correlation)

##### 4. Visualization
p <- ggplot(dt, aes(x = scale_rank_1, y = scale_rank_2, color = overlap_TF)) +
     geom_point(size = 2, alpha = 0.8) +
     # annotate("text", x = 0.8, y = 0.13, label = paste0("Correlation coefficient: \n", round(correlation, 3)), size = 5, fontface = "bold", color = "black") + 
     scale_color_manual(values = c("grey70", "#E41A1C")) + 
     theme_bw() + 
     theme(panel.grid.major = element_line(color = "grey85", linetype = "dotted"), panel.grid.minor = element_blank(), legend.position = "bottom") + 
     labs(x = "TF scaled rank (Geneset)", y = "TF scaled rank (Regionset)", color = "TF Overlap", title = paste0("Correlation coefficient: ", round(correlation, 3)))
ggsave("~/results/spRNA_spATAC_upstream_corr.pdf", p, width = 3, height = 3.5)

