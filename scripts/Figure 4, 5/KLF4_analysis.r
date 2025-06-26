
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)

spRNA_bart_result <- read.csv("~/results/spRNA/BART_results.csv") %>% mutate(rank = X) %>% filter(rank_avg_z_p_a_irwinhall_pvalue < 0.01) #35
spRNA_bart_result$TF <- str_to_title(tolower(spRNA_bart_result$TF))

subset <- readRDS("~/data/spRNA_subset.rds")
DefaultAssay(subset) <- "Spatial"

##### 1. Check expression of top-ranked TF coding gene
for (i in spRNA_bart_result$TF) {
     if (!i %in% rownames(subset)) {
        message(paste("Skipping gene:", i, " - not found"))
        next
     }

     expr_values <- FetchData(subset, vars = i, slot = "data")[[1]]
     upper_lim <- quantile(expr_values, 0.99, na.rm = TRUE)
     upper_lim <- 0.05

     p <- SpatialFeaturePlot(subset, features = i, pt.size.factor = 1.5, slot = "data") + 
            scale_fill_gradientn(colors = c("#FFF6F6", "#FFCCCC", "#FF9999", "#FF3333", "#FF0000"), limits = c(0, upper_lim), oob = scales::squish) + 
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.5, "cm"))
     ggsave(paste0("~/results/spRNA/top_TF_exp/", i, "_exp.pdf"), p, width = 5, height = 4)
     
     message("done: ", i)
}


##### 2. Check whether the coding gene belong to SVG or traj-gene
expression_matrix <- subset@assays$Spatial@counts
in_TF_name <- spRNA_bart_result[spRNA_bart_result$TF %in% rownames(expression_matrix), ] #35
not_in_TF_name <- spRNA_bart_result$TF[!spRNA_bart_result$TF %in% rownames(expression_matrix)] #13

traj_DEG_all <- read.csv("~/results/spRNA/traj_all_16159_genes.csv") 
traj_DEG <- traj_DEG_all %>% filter(adjusted_pvals < 0.1 & abs(correlation_rho) > 0.0921616)
traj_TF <- traj_DEG[traj_DEG$significant_features %in% spRNA_bart_result$TF, ] #15/35 = 42.8%

moran_DEG_all <- read.csv("~/results/spRNA/sp_all_10379_genes.csv") 
moran_DEG <- moran_DEG_all %>% filter(adjusted_pvals < 0.1 & Deviation_from_expectation > 0.07973)
sp_TF <- moran_DEG[moran_DEG$significant_features %in% spRNA_bart_result$TF, ] #15/35 = 42.8%

overlap_TF <- intersect(traj_TF$significant_features, sp_TF$significant_features) #9/35 = 25.7%
traj_DEG[traj_DEG$significant_features %in% overlap_TF, ]
moran_DEG[moran_DEG$significant_features %in% overlap_TF, ]

test1 <- traj_TF[!traj_TF$significant_features %in% overlap_TF, ]
test2 <- sp_TF[!sp_TF$significant_features %in% overlap_TF, ]
moran_DEG_all[moran_DEG_all$significant_features %in% test1$significant_features, ]
traj_DEG_all[traj_DEG_all$significant_features %in% test2$significant_features, ]

in_TF <- in_TF_name[in_TF_name$TF %in% overlap_TF, ] #9
in_TF <- in_TF %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)

not_in_TF <- in_TF_name[!in_TF_name$TF %in% c(traj_TF$significant_features, sp_TF$significant_features), ] #14
not_in_TF <- not_in_TF %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue, rank)




##### 3. Fisher exact test to test the motif spatial distribution at UDHS-level
library(GenomicRanges)

# UDHS
mm10_UDHS <- read.table("~/BART_library/mm10_library/bart2_mm10_UDHS.bed")
mm10_UDHS_gr <- GRanges(seqnames = mm10_UDHS$V1, ranges = IRanges(start = mm10_UDHS$V2 + 1, end = mm10_UDHS$V3)) #1529448

# spatially variable peaks on UDHS
moran_DAR <- read.csv("~/results/spATAC/sp_regionset.csv")
sp_peaks <- moran_DAR$significant_features #807
peak_info <- do.call(rbind, strsplit(sp_peaks, "-"))
output_df <- data.frame(V1 = peak_info[,1], V2 = as.integer(peak_info[,2]), V3 = as.integer(peak_info[,3]))
rownames(output_df) <- paste0(output_df$V1, "-", output_df$V2, "-", output_df$V3)
sp_peaks_gr <- GRanges(seqnames = output_df$V1, ranges = IRanges(start = output_df$V2, end = output_df$V3))
sp_peaks_on_UDHS <- findOverlaps(sp_peaks_gr, mm10_UDHS_gr, minoverlap = 50)
sp_overlap_indices <- subjectHits(sp_peaks_on_UDHS)
sp_peaks_UDHS <- mm10_UDHS_gr[sp_overlap_indices] #2102

# spatially non-variable peaks on UDHS
non_sp_overlap_indices <- setdiff(seq_along(mm10_UDHS_gr), sp_overlap_indices)
non_sp_peaks_UDHS <- mm10_UDHS_gr[non_sp_overlap_indices] #1526910

# KLF4 motif sites
bed_path <- paste0("~/motif_sites/mappable_results/KLF4_map.bed")
bed_files <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(bed_files) <- c("chr", "start", "end", "seq", "pval", "strand")
KLF4_motif_sites <- GRanges(seqnames = bed_files$chr, ranges = IRanges(start = bed_files$start + 1, end = bed_files$end)) # 298391
KLF4_motif_sites_on_UDHS <- findOverlaps(mm10_UDHS_gr, KLF4_motif_sites, minoverlap = 15)
KLF4_motif_sites_overlap_indices <- queryHits(KLF4_motif_sites_on_UDHS)
KLF4_motif_sites_UDHS <- mm10_UDHS_gr[KLF4_motif_sites_overlap_indices] # 32851


# Fisher's exact test
# Check overlap with KLF4 motif sites
sp_with_KLF4 <- findOverlaps(KLF4_motif_sites_UDHS, sp_peaks_UDHS, type = "within")
nonsp_with_KLF4 <- findOverlaps(KLF4_motif_sites_UDHS, non_sp_peaks_UDHS, type = "within") 

# Create counts for the 2x2 contingency table
a <- length(sp_with_KLF4)                                         # Spatial peaks with KLF4 motif
b <- length(sp_peaks_UDHS) - length(sp_with_KLF4)                 # Spatial peaks without KLF4 motif
c <- length(nonsp_with_KLF4)                                      # Non-spatial peaks with KLF4 motif
d <- length(non_sp_peaks_UDHS) - length(nonsp_with_KLF4)          # Non-spatial peaks without KLF4 motif

contingency_table <- matrix(c(a, b, c, d), nrow = 2, dimnames = list(Spatial = c("KLF4+", "KLF4-"), PeakType = c("Spatial", "Non-Spatial")))
fisher_result <- fisher.test(contingency_table, alternative = "greater")
print(fisher_result)

# Visualization
library(ggplot2)

plot_df <- data.frame(PeakType = rep(c("Spatial", "Non-Spatial"), each = 2),
                      KLF4 = rep(c("KLF4+", "KLF4-"), 2),
                      Count = c(a, b, c, d), 
                      Percentage = c(a/length(sp_peaks_UDHS)*100, b/length(sp_peaks_UDHS)*100, c/length(non_sp_peaks_UDHS)*100, d/length(non_sp_peaks_UDHS)*100))

p_val <- "< 2.2e-16"

p <- ggplot(plot_df, aes(x = PeakType, y = Percentage, fill = PeakType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  labs(title = paste0("Fisher's p-value = ", p_val), y = "Percentage", x = "Region Type") +
  scale_fill_manual(values = c("Spatial" = "#F47F72", "Non-Spatial" = "#F6B3AC")) +
  # scale_fill_manual(values = c("KLF4+" = "#F47F72", "KLF4-" = "#F6B3AC")) +
  # labs(title = paste0("Fisher's p-value = ", p_val), y = "Percentage", x = "Region Type", fill = "KLF4 Motif") +
  theme_bw() +
  ylim(0, 15) +
  theme(plot.title = element_text(size = 12), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave("~/results/spRNA/KLF4_sp_UDHS.pdf", p, width = 4, height = 2)
