
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

BART_downstream_TF_all <- read.csv("~/BART_results_downstream.csv") %>% mutate(rank = X)
BART_downstream_TF_all$TF <- str_to_title(tolower(BART_downstream_TF_all$TF))
spatrack_downstream_TF_all <- read.csv("~/spatrack_results_downstream.csv")
spatrack_downstream_TF_all <- spatrack_downstream_TF_all[!duplicated(spatrack_downstream_TF_all$TF), ] %>% mutate(rank = seq_len(nrow(.)))
scripro_downstream_TF_all <- read.csv("~/scripro_results_downstream.csv") %>% arrange(desc(downstream_activity)) %>% mutate(rank = seq_len(nrow(.)))
scripro_downstream_TF_all$TF <- scripro_downstream_TF_all$X

tf_interest <- c("Hnf4g", "Hnf4a", "Gata6", "Hnf1b", "Maf", "Cdx2", "Mafb", "Gata4", "Hes1")

BART_downstream_TF_all <- scale_rank(BART_downstream_TF_all, tf_interest)
spatrack_downstream_TF_all <- scale_rank(spatrack_downstream_TF_all, tf_interest)
scripro_downstream_TF_all <- scale_rank(scripro_downstream_TF_all, tf_interest)

BART_downstream_TF_all$method <- "BART-spatial"
spatrack_downstream_TF_all$method <- "SpaTrack"
scripro_downstream_TF_all$method <- "SCRIPro"

combined_rank <- bind_rows(BART_downstream_TF_all %>% select(TF, scaled_rank, TF_type, method),
                           spatrack_downstream_TF_all %>% select(TF, scaled_rank, TF_type, method), 
                           scripro_downstream_TF_all  %>% select(TF, scaled_rank, TF_type, method))

p10 <- VS_boxplot(combined_rank)
ggsave("benchmarking_boxplot.png", p10, width = 6, height = 3.5, dpi = 600)



cutoff_up <- 25
cutoff_down <- 25

BART_upstream_TF <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X) %>% head(., cutoff_up) 
BART_upstream_TF$TF <- str_to_title(tolower(BART_upstream_TF$TF))
BART_downstream_TF <- read.csv("~/BART_results_downstream.csv") %>% mutate(rank = X) %>% head(., cutoff_down) 
BART_downstream_TF$TF <- str_to_title(tolower(BART_downstream_TF$TF))
BART_TF <- c(BART_upstream_TF$TF, BART_downstream_TF$TF)
spatrack_upstream_TF <- read.csv("~/spatrack_results_upstream.csv") 
spatrack_upstream_TF <- spatrack_upstream_TF[!duplicated(spatrack_upstream_TF$TF), ] %>% mutate(rank = seq_len(nrow(.))) %>% head(., cutoff_up) 
spatrack_downstream_TF <- read.csv("~/spatrack_results_downstream.csv")
spatrack_downstream_TF <- spatrack_downstream_TF[!duplicated(spatrack_downstream_TF$TF), ] %>% mutate(rank = seq_len(nrow(.))) %>% head(., cutoff_down) 
spatrack_TF <- c(spatrack_upstream_TF$TF, spatrack_downstream_TF$TF)
scripro_upstream_TF <- read.csv("~/scripro_results_upstream.csv") %>% arrange(desc(upstream_activity)) %>% mutate(rank = seq_len(nrow(.))) %>% filter(p.value < 0.001) %>% head(., cutoff_up)
scripro_upstream_TF$TF <- scripro_upstream_TF$X
scripro_downstream_TF <- read.csv("~/scripro_results_downstream.csv") %>% arrange(desc(downstream_activity)) %>% mutate(rank = seq_len(nrow(.))) %>% filter(p.value < 0.001) %>% head(., cutoff_down)
scripro_downstream_TF$TF <- scripro_downstream_TF$X
scripro_TF <- c(scripro_upstream_TF$TF, scripro_downstream_TF$TF)

BART_TF_sig <- length(unique(BART_TF[BART_TF %in% tf_interest]))
spatrack_TF_sig <- length(unique(spatrack_TF[spatrack_TF %in% tf_interest]))
scripro_TF_sig <- length(unique(scripro_TF[scripro_TF %in% tf_interest]))

dt <- data.frame(method = c("BART-spatial", "SpaTrack", "SCRIPro"),
                 n_TFs  = c(BART_TF_sig, spatrack_TF_sig, scripro_TF_sig), 
                 in_TF = c("HNF4G\nHNF4A\nGATA6\nHNF1B", "MAF", "HNF4G\nHNF4A"))

p11 <- VS_barplot(dt)
ggsave("benchmarking_barplot.png", p11, width = 4, height = 3.5, dpi = 600)
