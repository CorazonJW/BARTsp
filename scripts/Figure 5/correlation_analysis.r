
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(ggrepel)

source("~/manuscript_figures.r")

##### 1. Load prediction results
mm_up_traj_only <- read.csv("~/BART_results_upstream_traj_only.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
mm_up_all <- read.csv("~/BART_results_upstream.csv") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)
mm_up_traj_only$scaled_rank <- round((max(mm_up_traj_only$rank) - mm_up_traj_only$rank) / (max(mm_up_traj_only$rank) - min(mm_up_traj_only$rank)), 3)
mm_up_all$scaled_rank <- round((max(mm_up_all$rank) - mm_up_all$rank) / (max(mm_up_all$rank) - min(mm_up_all$rank)), 3)

sig_dt <- merge(mm_up_traj_only, mm_up_all, by = "TF") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue.x, scaled_rank.x, rank_avg_z_p_a_irwinhall_pvalue.y, scaled_rank.y)
colnames(sig_dt) <- c("TF", "pseudo_time_only_pval", "pseudo-time_only_scaled_rank", "all_SVFs_pval", "all_SVFs_scaled_rank")
dt <- sig_dt

colnames(dt) <- c("TF", "pval1", "scale_rank_1", "pval2", "scale_rank_2")



##### 2. Test correlation
dt <- dt %>%
  mutate(color = case_when(pval1 < 0.05 & pval2 < 0.05 ~ "Shared TRs",   # shared
                           pval1 >= 0.05 & pval2 >= 0.05 ~ "Not significant",   # not sig
                           pval1 < 0.05 & pval2 >= 0.05 ~ "TVFs only unique",   # unique to pseudotime
                           pval1 >= 0.05 & pval2 < 0.05 ~ "Overlapped genes unique"))     # unique to all_SVFs

correlation_res <- cor.test(dt$scale_rank_1, dt$scale_rank_2, method = "spearman")
correlation <- round(correlation_res$estimate, 3)
print(correlation)



##### 3. Visualization
p4 <- ggplot(dt, aes(x = scale_rank_1, y = scale_rank_2, color = color)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("Shared TRs" = "#FCD807", "Not significant" = "#B4B4B4", 
                                "TVFs only unique" = "#51B14E", "Overlapped genes unique" = "#387EB8")) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey85", linetype = "dotted"), 
        panel.grid.minor = element_blank(), legend.position = "bottom", 
        legend.text = element_text(size = 14, family = "Arial"), 
        legend.title = element_blank(), 
        legend.direction = "vertical", 
        legend.key.width = unit(0.1, "cm"), 
        legend.spacing.x = unit(0.1, "cm"),
        legend.justification = c(0.5, 0),
        legend.box.spacing = unit(3, "pt"),
        legend.margin = margin(2,2,2,2),
        axis.text = element_text(size = 14, family = "Arial"), 
        axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
  labs(x = "TR scaled rank (TVFs only)", 
       y = "TR scaled rank (Overlapped features)", 
       color = "TR Category") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
ggsave("correlation.png", p4, width = 5, height = 5.5, dpi = 600)

