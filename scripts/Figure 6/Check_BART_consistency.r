
library(dplyr)

source("~/manuscript_figures.r")

spRNA <- read.csv("~/BART_results_upstream_RNA.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
spATAC <- read.csv("~/BART_results_upstream_ATAC.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
spATAC_rep1 <- read.csv("~/rep1/BART_results_upstream_ATAC.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
spATAC_rep2 <- read.csv("~/rep2/BART_results_upstream_ATAC.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)

p1 <- correlation(spRNA, spATAC, "TR scaled ranks (Geneset)", "TR scaled ranks (Regionset)")
p2 <- correlation(spATAC, spATAC_rep1, "TR scaled ranks (Regionset)", "TR scaled ranks (spATAC rep1)")
p3 <- correlation(spATAC, spATAC_rep2, "TR scaled ranks (Regionset)", "TR scaled ranks (spATAC rep2)")
p4 <- correlation(spATAC_rep1, spATAC_rep2, "TR scaled ranks (spATAC rep1)", "TR scaled ranks (spATAC rep2)")

ggsave("spRNA_spATAC.png", p1, height = 4, width = 4)
ggsave("spATAC_spATAC1.png", p2, height = 4, width = 4)
ggsave("spATAC_spATAC2.png", p3, height = 4, width = 4)
ggsave("spATAC1_spATAC2.png", p4, height = 4, width = 4)
