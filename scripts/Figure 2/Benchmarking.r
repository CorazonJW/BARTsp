
source("~/manuscript_figures.r")

library(BARTsp)
library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)
library(sciColors)
library(forcats)
library(stringr)
library(gridExtra)

palette <- c("BART-spatial" = "#4C78A8", "SpaTrack" = "#F58518",  "SCRIPro" = "#54A24B")

tf_interest <- c("GRHL3", "TCF4", "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", "FOS", "JUN", "JUND")

BART_upstream_TF_all <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X)
spatrack_upstream_TF_all <- read.csv("~/spatrack_results_upstream.csv") 
spatrack_upstream_TF_all <- spatrack_upstream_TF_all[!duplicated(spatrack_upstream_TF_all$TF), ] %>% mutate(rank = seq_len(nrow(.))) 
scripro_upstream_TF_all <- read.csv("~/scripro_results_upstream.csv") %>% arrange(desc(core_activity)) %>% mutate(rank = seq_len(nrow(.)))
scripro_upstream_TF_all$TF <- scripro_upstream_TF_all$X

BART_upstream_TF_all <- scale_rank(BART_upstream_TF_all, tf_interest)
spatrack_upstream_TF_all <- scale_rank(spatrack_upstream_TF_all, tf_interest)
scripro_upstream_TF_all <- scale_rank(scripro_upstream_TF_all, tf_interest)

BART_upstream_TF_all$method <- "BART-spatial"
spatrack_upstream_TF_all$method <- "SpaTrack"
scripro_upstream_TF_all$method <- "SCRIPro"

combined_rank <- bind_rows(BART_upstream_TF_all %>% select(TF, scaled_rank, TF_type, method),
                           spatrack_upstream_TF_all %>% select(TF, scaled_rank, TF_type, method), 
                           scripro_upstream_TF_all  %>% select(TF, scaled_rank, TF_type, method))

p10 <- VS_boxplot(combined_rank)
ggsave("benchmarking_boxplot.png", p10, width = 6, height = 3.5, dpi = 600)



cutoff <- 25

BART_upstream_TF <- BART_upstream_TF_all %>% head(., cutoff) 
spatrack_upstream_TF <- spatrack_upstream_TF_all[!duplicated(spatrack_upstream_TF_all$TF), ] %>% mutate(rank = seq_len(nrow(.))) %>% head(., cutoff) 
scripro_upstream_TF <- scripro_upstream_TF_all %>% arrange(desc(core_activity)) %>% mutate(rank = seq_len(nrow(.))) %>% filter(p.value < 0.01) %>% head(., cutoff)
scripro_upstream_TF$TF <- scripro_upstream_TF_all$X

BART_TF_sig <- length(unique(BART_upstream_TF[BART_upstream_TF$TF %in% tf_interest, ]$TF))
spatrack_TF_sig <- length(unique(spatrack_upstream_TF[spatrack_upstream_TF$TF %in% tf_interest, ]$TF))
scripro_TF_sig <- length(unique(scripro_upstream_TF[scripro_upstream_TF$TF %in% tf_interest, ]$TF))

dt <- data.frame(method = c("BART-spatial", "SpaTrack", "SCRIPro"),
                 n_TFs  = c(BART_TF_sig, spatrack_TF_sig, scripro_TF_sig),
                 in_TF = c("KLF4\nGRHL2\nTP63\nFOS\nSNAI2\nJUN", "GRHL3", "STAT3\nKLF4\nFOS\nTP63\nJUND"))

p11 <- VS_barplot(dt)
ggsave("benchmarking_barplot.png", p11, width = 4, height = 3.5, dpi = 600)

