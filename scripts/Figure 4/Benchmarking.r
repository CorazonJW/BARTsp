
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)
library(stringr)
library(ggsignif)
library(forcats)
library(stringr)

source("~/manuscript_figures.r")

cutoff_up <- 25
cutoff_down <- 25

BART_upstream_TF <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01) %>% head(., cutoff_up) 
BART_upstream_TF$TF <- str_to_title(tolower(BART_upstream_TF$TF))
BART_downstream_TF <- read.csv("~/BART_results_downstream.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01) %>% head(., cutoff_down)
BART_downstream_TF$TF <- str_to_title(tolower(BART_downstream_TF$TF))

spatrack_upstream_TF <- read.csv("~/spatrack_results_upstream.csv")
spatrack_upstream_TF <- spatrack_upstream_TF[!duplicated(spatrack_upstream_TF$TF), ] %>% mutate(rank = seq_len(nrow(.))) %>% head(., cutoff_up) 
spatrack_downstream_TF <- read.csv("~/spatrack_results_downstream.csv")
spatrack_downstream_TF <- spatrack_downstream_TF[!duplicated(spatrack_downstream_TF$TF), ] %>% mutate(rank = seq_len(nrow(.))) %>% head(., cutoff_down) 

scripro_upstream_TF <- read.csv("~/scripro_results_upstream.csv") %>% arrange(desc(upstream_activity)) %>% mutate(rank = seq_len(nrow(.))) %>% filter(p.value < 0.01) %>% head(., cutoff_up)
scripro_upstream_TF$TF <- scripro_upstream_TF$X
scripro_downstream_TF <- read.csv("~/scripro_results_downstream.csv") %>% arrange(desc(downstream_activity)) %>% mutate(rank = seq_len(nrow(.))) %>% filter(p.value < 0.01) %>% head(., cutoff_down)
scripro_downstream_TF$TF <- scripro_downstream_TF$X


palette <- c("BART-spatial\n(upstream)" = "#4C78A8", "BART-spatial\n(downstream)" = "#3167A5", 
             "SpaTrack\n(upstream)" = "#F0A253", "SpaTrack\n(downstream)" = "#F58518", 
             "SCRIPro\n(upstream)" = "#54A24B", "SCRIPro\n(downstream)" = "#2B9B1E")


tf_interest_up <-  c("Pax6", "Sox9", "Fezf2", "Klf4")
tf_interest_down <-  c("Neurod2")

BART_TF_sig_up <- length(unique(BART_upstream_TF[BART_upstream_TF$TF %in% tf_interest_up, ]$TF))
spatrack_TF_sig_up <- length(unique(spatrack_upstream_TF[spatrack_upstream_TF$TF %in% tf_interest_up, ]$TF))
scripro_TF_sig_up <- length(unique(scripro_upstream_TF[scripro_upstream_TF$TF %in% tf_interest_up, ]$TF))
BART_TF_sig_down <- length(unique(BART_downstream_TF[BART_downstream_TF$TF %in% tf_interest_down, ]$TF))
spatrack_TF_sig_down <- length(unique(spatrack_downstream_TF[spatrack_downstream_TF$TF %in% tf_interest_down, ]$TF))
scripro_TF_sig_down <- length(unique(scripro_downstream_TF[scripro_downstream_TF$TF %in% tf_interest_down, ]$TF))


dt <- data.frame(method = c("BART-spatial\n(upstream)", "SpaTrack\n(upstream)", "SCRIPro\n(upstream)", 
                            "BART-spatial\n(downstream)",  "SpaTrack\n(downstream)", "SCRIPro\n(downstream)"),
                 n_TFs  = c(BART_TF_sig_up, spatrack_TF_sig_up, scripro_TF_sig_up, 
                            BART_TF_sig_down, spatrack_TF_sig_down, scripro_TF_sig_down), 
                 in_TF = c("SOX9\nFEZF2\nKLF4", "", "SOX9", 
                           "NEUROD2", "", ""))

dt$method <- factor(dt$method, levels = c("BART-spatial\n(upstream)", "SpaTrack\n(upstream)", "SCRIPro\n(upstream)", 
                                          "BART-spatial\n(downstream)",  "SpaTrack\n(downstream)", "SCRIPro\n(downstream)"))

p5 <- VS_barplot(dt) + geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", size = 0.6)
ggsave("benchmarking_barplot.png", p5, width = 8, height = 3.5, dpi = 600)








BART_upstream_TF_all <- read.csv("~/BART_results_upstream.csv") %>% mutate(rank = X)
BART_upstream_TF_all$TF <- str_to_title(tolower(BART_upstream_TF_all$TF))
BART_downstream_TF_all <- read.csv("~/BART_results_downstream.csv") %>% mutate(rank = X) 
BART_downstream_TF_all$TF <- str_to_title(tolower(BART_downstream_TF_all$TF))

spatrack_upstream_TF_all <- read.csv("~/spatrack_results_upstream.csv")
spatrack_upstream_TF_all <- spatrack_upstream_TF_all[!duplicated(spatrack_upstream_TF_all$TF), ] %>% mutate(rank = seq_len(nrow(.)))
spatrack_downstream_TF_all <- read.csv("~/spatrack_results_downstream.csv") %>% mutate(rank = seq_len(nrow(.)))
spatrack_downstream_TF_all <- spatrack_downstream_TF_all[!duplicated(spatrack_downstream_TF_all$TF), ] %>% mutate(rank = seq_len(nrow(.)))

scripro_upstream_TF_all <- read.csv("~/scripro_results_upstream.csv") %>% arrange(desc(upstream_activity)) %>% mutate(rank = seq_len(nrow(.))) 
scripro_upstream_TF_all$TF <- scripro_upstream_TF_all$X
scripro_downstream_TF_all <- read.csv("~scripro_results_downstream.csv") %>% arrange(desc(downstream_activity)) %>% mutate(rank = seq_len(nrow(.)))
scripro_downstream_TF_all$TF <- scripro_downstream_TF_all$X

tf_interest_up <-  c("Pax6", "Sox9", "Fezf2", "Klf4")

BART_upstream_TF_all <- scale_rank(BART_upstream_TF_all, tf_interest_up)
spatrack_upstream_TF_all <- scale_rank(spatrack_upstream_TF_all, tf_interest_up)
scripro_upstream_TF_all <- scale_rank(scripro_upstream_TF_all, tf_interest_up)

BART_upstream_TF_all$method <- "BART-spatial"
spatrack_upstream_TF_all$method <- "SpaTrack"
scripro_upstream_TF_all$method <- "SCRIPro"

combined_rank <- bind_rows(
  BART_upstream_TF_all     %>% select(TF, scaled_rank, TF_type, method),
  spatrack_upstream_TF_all %>% select(TF, scaled_rank, TF_type, method), 
  scripro_upstream_TF_all  %>% select(TF, scaled_rank, TF_type, method)
)

p6 <- VS_boxplot(combined_rank)
ggsave("benchmarking_boxplot.png", p2, width = 6, height = 3.5, dpi = 600)