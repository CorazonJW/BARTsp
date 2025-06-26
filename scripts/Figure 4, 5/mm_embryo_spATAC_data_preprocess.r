
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ArchR)
library(grid)

setwd("~/spATAC_output/")
data_dir <- "~/data/"
output_dir <- "~/spATAC_output/"

##### 1. Pre-process
# Load data
threads = 8
addArchRThreads(threads = threads)
addArchRGenome("mm10")
inputFiles <- paste0(data_dir, "peak_expression_matrix.tsv.gz")
sampleNames <- 'ME13'
# Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000), 
  force = TRUE
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE, 
  force = TRUE
)
# Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data_dir, "ME13_50um_spatial"), filter.matrix = filter.matrix)
name <- paste0("ME13#", rownames(image@coordinates), "-1")
rownames(image@coordinates) <- name
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue
# Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)
proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)
proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
proj_in_tissue <- addImputeWeights(proj_in_tissue)


##### 2. Add MOCA predicted id to ArchR object
E13_sp <- readRDS("~/spRNA_output/spRNA_preprocess.rds")

dt <- data.frame(E13_sp$predicted.id)
rownames(dt) <- paste0("ME13#", rownames(dt), "-1")
dt <- dt[match(proj_in_tissue@cellColData@rownames, rownames(dt)), , drop = FALSE]

proj_in_tissue$predicted.id <- dt$E13_sp.predicted.id


##### 3. Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "predicted.id")

pathToMacs2 <- findMacs2()

proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "predicted.id", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

saveRDS(proj_in_tissue, paste0(output_dir, "spatac_preprocess.RDS"))
proj_in_tissue <- readRDS(paste0(output_dir, "spatac_preprocess.RDS"))



###### 4. Motif analysis to add motif information to ArchR object
# ChromVAR Deviatons Enrichment
if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
    proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif")
}

proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveRDS(proj_in_tissue, paste0(output_dir, "spatac_preprocess.RDS"))
proj_in_tissue <- readRDS(paste0(output_dir, "spatac_preprocess.RDS"))


# Subset project to only contain radial glia and postmitotic premature neurons
E13_spRNA <- readRDS("~/spRNA_output/spRNA_preprocess.rds")
Idents(E13_spRNA) <- E13_spRNA$predicted.id
subset <- subset(E13_spRNA, idents = c("Radial glia", "Postmitotic premature neurons"))
cell_subset <- paste0("ME13#", rownames(subset@meta.data), "-1")
ct_label <- as.character(subset$predicted.id)

subset_ATAC <- subsetArchRProject(
  ArchRProj = proj_in_tissue,
  cells = cell_subset,
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = subset_ATAC, 
  useMatrix = "MotifMatrix", 
  groupBy = "predicted.id",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)

saveRDS(subset_ATAC, paste0(output_dir, "spatac_subset.RDS"))
subset_ATAC <- readRDS(paste0(output_dir, "spatac_subset.RDS"))

# ArchR Motif enrichment
subset_ATAC <- readRDS(paste0(output_dir, "spatac_subset.RDS"))

markerTest <- getMarkerFeatures(
  ArchRProj = subset_ATAC, 
  useMatrix = "PeakMatrix",
  groupBy = "predicted.id",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = subset_ATAC,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 3,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(paste0(output_dir, "Motif.pdf"), ggDo, width = 5, height = 5)



##### 5. Construct trajectory
subset_ATAC <- readRDS(paste0(output_dir, "spatac_subset.RDS"))
trajectory <- c("Radial glia", "Postmitotic premature neurons")

proj_traj <- addTrajectory(
  ArchRProj = subset_ATAC, 
  name = "trajectory", 
  groupBy = "predicted.id",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
head(proj_traj$trajectory[!is.na(proj_traj$trajectory)])

saveRDS(proj_traj, paste0(output_dir, "spatac_trajectory.RDS"))
proj_traj <- readRDS(paste0(output_dir, "spatac_trajectory.RDS"))

# Plot Motif
plotVarDev <- getVarDeviations(proj_traj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_traj, addDOC = FALSE)

# Visualize spATAC trajectory 
trajMM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajPM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)

# spatac_gene <- sub("chr[0-9XY]+:", "", trajGSM@NAMES)
# spatac_motif <- gsub("^(deviations:|z:)|_[0-9]+$", "", trajMM@NAMES)
# spatac_peak <- gsub("chr(\\d+):(\\d+)_(\\d+)", "chr1-\\2-\\3", trajPM@NAMES)

pdf(paste0(output_dir, "pseudotime_motif.pdf"), width = 10, height = 8)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
p1
p2
p3
dev.off()

# Visualize spATAC trajectory 
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")

dt <- data.frame(cell = proj_traj@cellColData@rownames, sp_ATAC_traj = proj_traj$trajectory)
dt$cell <- gsub("^ME13#|-1$", "", dt$cell)
dt <- dt[match(WhichCells(subset), dt$cell), ]

subset <- AddMetaData(object = subset, metadata = dt$sp_ATAC_traj, col.name = "spATAC_traj")

DefaultAssay(subset) <- "ATAC"
subset <- subset[, !is.na(subset$spATAC_traj)]

p1 <- SpatialFeaturePlot(subset, features = "pseudotime", alpha = c(0.9, 0.9), pt.size.factor = 1.5) + guides(fill = guide_colorbar(title = "Pseudotime (spRNA)")) +
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.5, "cm"))
p2 <- SpatialFeaturePlot(subset, features = "spATAC_traj", alpha = c(0.9, 0.9), pt.size.factor = 1.5) + guides(fill = guide_colorbar(title = "Pseudotime (spATAC)")) +
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.5, "cm"))
p <- p1+p2
ggsave(paste0(output_dir, "spATAC_traj.pdf"), p, width = 10, height = 4)

saveRDS(subset, "/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/E13/spRNA_output/spRNA_subset.rds")










##### 6. Trajectory analysis: Check peaks whose activity correlate with trajectory
# Method 1: Use peak matrix and activity change
proj_traj <- readRDS(paste0(output_dir, "spatac_trajectory.RDS"))
trajPM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)

z_mat <- assays(trajPM)$smoothMat  
rownames(z_mat) <- sub(":", "-", rownames(z_mat))
rownames(z_mat) <- sub("_", "-", rownames(z_mat))

# Compute row means for each segment
intervals <- seq(1, 100, by = 10) 
mean_values <- lapply(intervals, function(start_col) {
  cols <- start_col:(start_col + 9)  # Select 10-column segment
  rowMeans(z_mat[, cols, drop = FALSE])  # Compute row means
})

mean_values_matrix <- do.call(cbind, mean_values)
colnames(mean_values_matrix) <- paste0("mean_", intervals)

# Identify rows where any two segments have a significant difference (≥ 0.8)
significant_change <- apply(mean_values_matrix, 1, function(row_means) {
  for (i in 1:(length(row_means) - 1)) {
    for (j in (i + 1):length(row_means)) {
      if (abs(row_means[i] - row_means[j]) >= 0.8) {
        return(TRUE)
      }
    }
  }
  return(FALSE) 
})

z_mat_filtered <- z_mat[significant_change, ]
str(z_mat_filtered)
peaks <- rownames(z_mat_filtered)
write.table(peaks, paste0(output_dir, "significant_peaks_trajectory_expression.txt"))


# Method 2: Use correlation test
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")

expression_matrix <- subset@assays$peaks@counts
pseudotime_values <- subset$spATAC_traj

cor_results <- apply(expression_matrix, 1, function(expr) {
    cor.test(expr, pseudotime_values, method = "spearman")$p.value
})

adjusted_pvals <- p.adjust(cor_results, method = "fdr")
significant_peak <- names(adjusted_pvals[adjusted_pvals < 0.05]) 
significant_peak <- significant_peak[!is.na(significant_peak)] # 1095 

write.table(significant_peak, paste0(output_dir, "significant_peaks_trajectory_correlation.txt"))


# Check the overlap
library(GenomicRanges)

expression_peaks <- read.table(paste0(output_dir, "significant_peaks_trajectory_expression.txt")) #675
expression_peaks <- expression_peaks$x
correlation_peaks <- read.table(paste0(output_dir, "significant_peaks_trajectory_correlation.txt")) #1095
correlation_peaks <- correlation_peaks$x

convert_to_GRanges <- function(peak_list) {
  split_peaks <- do.call(rbind, strsplit(peak_list, "-"))
  GRanges(seqnames = split_peaks[, 1], 
          ranges = IRanges(start = as.numeric(split_peaks[, 2]), 
                           end = as.numeric(split_peaks[, 3])))
}

gr_exp_peak <- convert_to_GRanges(expression_peaks)
gr_corr_peak <- convert_to_GRanges(correlation_peaks)
overlaps <- findOverlaps(gr_exp_peak, gr_corr_peak, maxgap = 200)

overlapping_peaks_1 <- gr_exp_peak[queryHits(overlaps)]
overlapping_peaks_2 <- gr_corr_peak[subjectHits(overlaps)]
overlap_peaks <- data.frame(Overlap_Peaks_1 = as.character(overlapping_peaks_1), Overlap_Peaks_2 = as.character(overlapping_peaks_2))
write.table(overlap_peaks, paste0(output_dir, "significant_peaks_trajectory_overlap.txt"))


# Method 3: Use variance test
proj_traj <- readRDS(paste0(output_dir, "spatac_trajectory.RDS"))
trajPM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)

z_mat <- assays(trajPM)$smoothMat  
rownames(z_mat) <- sub(":", "-", rownames(z_mat))
rownames(z_mat) <- sub("_", "-", rownames(z_mat))
colnames(z_mat) <- as.character(c(1:100))

variance <- apply(z_mat, 1, var)

# visualize variance
library(ggplot2)
library(patchwork)

variance_df <- data.frame(variance = as.numeric(variance))
p1 <- ggplot(variance_df, aes(x = variance)) +
  geom_density(fill = "grey", alpha = 0.8) +
  labs(x = "Variance", y = "Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2 <- ggplot(variance_df, aes(y = variance)) +
  geom_boxplot(alpha = 0.8, size = 0.3, outlier.size = 0.5) +
  labs(y = "Variance") +
  geom_hline(yintercept = quantile(variance, 0.75) + 1.5 * (quantile(variance, 0.75) - quantile(variance, 0.25)), linetype = "dashed", color = "red", size = 0.3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p1 / p2 
ggsave("~/spATAC_output/traj_sig_peak_var_test.pdf", p, width = 2, height = 4)

# use standard deviation cutoff (1470 peaks)
cutoff_pos <- mean(variance) + 2 * sd(variance)  # Mean + 2 SD
cutoff_neg <- mean(variance) - 2 * sd(variance)  # Mean - 2 SD
significant_peak_sd <- variance[variance > cutoff_pos | variance < cutoff_neg]

# use IQR cutoff to select significant peaks (more strict, 1147 peaks)
IQR <- quantile(variance, 0.75) - quantile(variance, 0.25)
cutoff_pos <- quantile(variance, 0.75) + 1.5 * IQR # Q3 + 1.5*IQR
cutoff_neg <- quantile(variance, 0.25) - 1.5 * IQR # Q1 - 1.5*IQR
significant_peak_IQR <- variance[variance > cutoff_pos | variance < cutoff_neg]

# use permutation method to set cutoff (most strict, 833 peaks)
set.seed(2025)
permuted_variances <- replicate(10, var(apply(z_mat, 1, sample)))
cutoff_permut <- quantile(permuted_variances, 0.95) 
significant_peak_permut <- variance[variance > cutoff_permut]

significant_peak <- names(significant_peak_IQR)
write.table(significant_peak, paste0(output_dir, "significant_peaks_trajectory_variance.txt"))


# Method 4: use generalized additive model
# 1. fits a smooth curve to model peak accessibility changes along pseudotime
# 2. uses non-parametric splines to detect peaks with dynamic trends
proj_traj <- readRDS(paste0(output_dir, "spatac_trajectory.RDS"))
trajPM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
z_mat <- assays(trajPM)$smoothMat  
rownames(z_mat) <- sub(":", "-", rownames(z_mat))
rownames(z_mat) <- sub("_", "-", rownames(z_mat))
colnames(z_mat) <- as.character(c(1:100))

write.table(z_mat, "~/spATAC_output/z_mat.txt")
z_mat <- read.table("~/spATAC_output/z_mat.txt")
z_mat <- as.matrix(z_mat)

library(mgcv)
gam_model <- lapply(seq_len(nrow(z_mat)), function(i) {
  df <- data.frame(pseudotime = c(1:100), accessibility = as.numeric(z_mat[i, ]))
  gam(accessibility ~ s(pseudotime, bs="cs"), data = df)
})

# extract pvalues for the smooth term
p_value <- vector("list", nrow(z_mat))

for (i in 1:nrow(z_mat)) {
  p_value[[i]] <- tryCatch({
    gam_summary <- summary(gam_model[[i]])
    gam_summary$s.table["s(pseudotime)", "p-value"]
  }, error = function(e) {
    message("Skipping row ", i, ": ", e$message)
    NA
  })
}

p_values <- unlist(p_value)
df_pvalues <- data.frame(peak_name = rownames(z_mat), p_value = p_values)
df_pvalues <- na.omit(df_pvalues)
df_pvalues$adj_p_value <- p.adjust(df_pvalues$p_value, method = "fdr")
df_pvalues <- df_pvalues[order(df_pvalues$adj_p_value), ]

significant_peak <- df_pvalues %>% filter(p_value < 0.01) # 31730 regions are signigicant


##### 7. Functional annotation of trajectory_reslated significant peaks for verification
library(GenomicRanges)
library(rGREAT)
expression_peaks <- read.table(paste0(output_dir, "significant_peaks_trajectory_expression.txt")) #675
expression_peaks <- expression_peaks$x
correlation_peaks <- read.table(paste0(output_dir, "significant_peaks_trajectory_correlation.txt")) #1095
correlation_peaks <- correlation_peaks$x
variance_peaks <- read.table(paste0(output_dir, "significant_peaks_trajectory_variance.txt")) # 1147
variance_peaks <- variance_peaks$x

# convert to gr object
traj_sig_peaks <- variance_peaks
gr_significant <- GRanges(
  seqnames = gsub("-.*", "", traj_sig_peaks), 
  ranges = IRanges(start = as.numeric(gsub(".*-(\\d+)-\\d+", "\\1", traj_sig_peaks)),
                   end = as.numeric(gsub(".*-\\d+-(\\d+)", "\\1", traj_sig_peaks)))
)

# run GREAT
res_BP <- great(gr_significant, "GO:BP", "TxDb.Mmusculus.UCSC.mm10.knownGene")
tb_BP <- getEnrichmentTable(res_BP)

BP_result <- tb_BP %>% select(description, genome_fraction, observed_region_hits, fold_enrichment, p_value, p_adjust) %>% filter(fold_enrichment >= 2 & p_adjust < 0.05)

BP_result_corr <- read.csv(paste0(output_dir, "traj_sig_peaks_corr_GO_result.csv"))
BP_result_exp <- read.csv(paste0(output_dir, "traj_sig_peaks_exp_GO_result.csv"))
BP_result_var <- read.csv(paste0(output_dir, "traj_sig_peaks_var_GO_result.csv"))

sig_pathway <- c("growth", "central nervous system development", "developmental growth", "brain development", "negative regulation of cell differentiation", 
                 "regulation of growth", "forebrain development", "regulation of developmental growth", "cell fate commitment", "neural precursor cell proliferation", 
                 "glial cell differentiation", "cerebral cortex development", "stem cell population maintenance", "cell proliferation in forebrain", 
                 "negative regulation of neuron differentiation", "cell fate specification", "neuron fate specification", "brain development", "glial cell development",
                 "positive regulation of developmental growth", "radial glial cell differentiation", "glial cell fate specification", "negative regulation of glial cell differentiation", 
                 "cell morphogenesis involved in neuron differentiation", "neuron projection development", "positive regulation of growth", "forebrain radial glial cell differentiation", 
                 "neuron fate commitment", "forebrain cell migration")

sig_pathway_result_corr <- BP_result_corr %>% filter(description %in% sig_pathway)
sig_pathway_result_exp <- BP_result_exp %>% filter(description %in% sig_pathway)
sig_pathway_result_var <- BP_result_var %>% filter(description %in% sig_pathway)
nrow(sig_pathway_result_corr)
nrow(sig_pathway_result_exp)
nrow(sig_pathway_result_var)

# visualization
library(ggplot2)
library(forcats)
p <- ggplot(BP_result_var_dt, aes(x = observed_region_hits, y = fct_reorder(description, observed_region_hits), color = fold_enrichment, size = observed_region_hits)) +
  geom_point() + 
  theme_bw() +
  labs(x = "Observed Region Hits", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", axis.text.y = element_text(size = 11)) +
  guides(color = guide_colorbar(title = "Fold enrichment"), size = guide_legend(title = "-log10(pval_adj)"))
ggsave(paste0(output_dir, "traj_sig_peaks_var_GO_result.pdf"), p, height = 4, width = 7)










###### 8. Spatailly variable peak analysis: find genes whose expression varies in radial glia and PPN region
# Method 0: Test with SPARK (not working)
library('SPARK')

spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_coor <- paste0(spatial_loc$x, "x", spatial_loc$y)
spatial_coor <- spatial_coor[order(match(spatial_coor, colnames(peak_count)))]

subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak_count <- subset@assays$peaks@counts
colnames(peak_count) <- spatial_coor

info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(peak_count),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(peak_count),split="x"),"[",2)))
rownames(info)  <- colnames(peak_count)

spark <- CreateSPARKObject(counts = peak_count, location = info[,1:2], percentage = 0.1, min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)

spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 4, verbose = F)
spark <- spark.test(spark, check_positive = T, verbose = F)
head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])

# sparkX <- sparkx(peak_count, location, numCores=1, option="mixture")


# Method 1: Test with Seurat Moran's I
library(Seurat)

# expression matrix
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts
peaks <- rownames(peak) 

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])

results <- RunMoransI(peak, S, verbose = TRUE)
sig_results <- test %>% filter(p.value < 0.05 & observed !=0)

write.table(sig_results, paste0(output_dir, "spatially_significant_peaks_seurat_moran's_I.txt"))



# Method 2: Test with Moran's I to find spatially differential expressed peaks
library(spdep)
library(Matrix)
library(expm)
library(stats)
library(dplyr)
library(parallel)

# Load data
# peak expression
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts
peaks <- rownames(peak) 

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])

# 1. Center and scale gene expression and spatial coordinate matrix
preprocess_data <- function(expression_matrix, spatial_coordinates) {
  y <- scale(expression_matrix, center = TRUE, scale = TRUE)
  S <- scale(spatial_coordinates, center = TRUE, scale = TRUE)

  return(list(y = y, S = S))
}

# 2. Compute Moran's I for spatial autocorrelation
compute_morans_I <- function(y, S) {
  coords <- as.data.frame(S)
  k <- min(5, nrow(coords) - 1)
  neighbors <- knearneigh(coords, k = k)  # Find nearest neighbors
  listw <- nb2listw(knn2nb(neighbors), style = "W")
  
  morans_test <- moran.test(y, listw)
  return(morans_test)
}

# 3. Process all peaks 
morana_I_result <- list()
for (i in peaks) {
    y <- as.matrix(peak[as.character(i), ])
    colnames(y) <- as.character(i)
    # Skip if all values in y are zero
    if (all(y == 0, na.rm = TRUE)) {
        message(paste("Skipping peak:", i, " - all values are zero"))
        next
    }
    # Skip if y contains NA values
    if (any(is.na(y))) {
        message(paste("Skipping peak:", i, " - contains NA values"))
        next
    }
    processed_data <- preprocess_data(y, S)
    # Skip if NA values remain after processing
    if (any(is.na(processed_data$y))) {
        message(paste("Skipping peak:", i, "- NA values in scaled expression"))
        next
    }
    result <- compute_morans_I(processed_data$y, processed_data$S)
    if (!is.null(result)) {
        morana_I_result[[i]] <- result
    }
}

# Adjust p-value with fdr
p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
  morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

saveRDS(morana_I_result, paste0(output_dir, "moran's_I_result.RDS"))

# Select significant peaks
significant_result <- morana_I_result[sapply(morana_I_result, function(x) x$estimate[[1]] != x$estimate[[2]])] #32293/32361
significant_result <- morana_I_result[sapply(morana_I_result, function(x) x$adjusted_p.value < 0.05)] #583

sp_sig_peaks <- names(significant_result)
write.table(sp_sig_peaks, paste0(output_dir, "spatially_significant_peaks.txt"))

positive_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] > x$estimate[[2]])] # all positive regions
negative_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] < x$estimate[[2]])]



# Method 3: Test with singleCellHaystack to find spatially differential expressed peaks (not very well)
set.seed(2025)

# peak expression
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])

# run singleCellHaystack
library(singleCellHaystack)
res <- haystack(S, peak)
result <- show_result_haystack(res.haystack = res) %>% filter(log.p.adj < -1.3)
write.table(result, paste0(output_dir, "spatially_significant_peaks_haystack.txt"))




# Method 4: Test with SPADE to find spatially differential expressed peaks (not very well)
# peak expression
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])
S <- as.data.frame(S)

# run SPADE
library(SPADE)
data_norm <- SPADE_norm(readcounts = as.matrix(peak), info = S)
Est <- SPADE_estimate(expr_data = data_norm, info = S)
Test_res <- SPADE_test(object = data_norm, location = S, para = Est)
result <- Test_res %>% filter(Adjust.Pvalue < 0.05)
write.table(result, paste0(output_dir, "spatially_significant_peaks_SPADE.txt"))




# Method 5: Test with SpaGene to find spatially differential expressed peaks
# peak expression
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])
S <- as.data.frame(S)

# run SpaGene
library(SpaGene)
mc_sv <- SpaGene(peak,S)
result <- mc_sv$spagene_res[order(mc_sv$spagene_res$adjp),] %>% filter(adjp < 0.05)
write.table(result, paste0(output_dir, "spatially_significant_peaks_SpaGene.txt"))





# Method 6: Test with MERINGUE to find spatially differential expressed peaks
# peak expression
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
peak <- subset@assays$peaks@counts

# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])
S <- as.data.frame(S)

# run MERINGUE
library(MERINGUE)
N <- getSpatialNeighbors(S)
results <- getSpatialPatterns(peak, N)
filter <- filterSpatialPatterns(mat = peak, I = results, w = N, alpha = 0.05, details = TRUE, minPercentCells = 0.05)
write.table(filter, paste0(output_dir, "spatially_significant_peaks_MERINGUE.txt"))

# check overlap of MERINGUE results and Moran's I test
MERINGUE_result <- read.table(paste0(output_dir, "spatially_significant_peaks_MERINGUE.txt"))
MERINGUE_peak <- rownames(MERINGUE_result) #158
moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks.txt"))
moran_peak <- moran_result$x #583

overlap_peak <- intersect(MERINGUE_peak, moran_peak) #98 (good overlap)

# check overlap of MERINGUE results, Seurat Moran's I test, and self-written Moran's I
seurat_moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks_seurat_moran's_I.txt"))
seurat_moran_peak <- rownames(seurat_moran_result) #2166
moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks.txt"))
moran_peak <- moran_result$x #583

overlap_peak <- intersect(seurat_moran_peak, moran_peak) #507 (good overlap)
overlap_peak <- intersect(MERINGUE_peak, seurat_moran_peak) #140 (good overlap)





##### 8. Look for peaks that are spatially significant and also correlate with trajectory
library(GenomicRanges)

# spatially important peaks
sp_sig_peaks <- read.table("~/spATAC_output/spatially_significant_peaks.txt")
sp_sig_peaks <- sp_sig_peaks$x # 583

# 1. overlap with peaks detected by correlation test
traj_sig_peaks_coor <- read.csv(paste0(output_dir, "significant_peaks_trajectory_correlation.txt"))
traj_sig_peaks_coor <- gsub("^\\d+\\s+", "", traj_sig_peaks_coor$x) # 1095

overlap_sig_peak <- intersect(sp_sig_peaks, traj_sig_peaks_coor) #101
write.table(overlap_sig_peak, paste0(output_dir, "overlap_peaks_sp_traj_corr.txt"))

# 2. overlap with peaks detected by expression
traj_sig_peaks_exp <- read.table(paste0(output_dir, "significant_peaks_trajectory_expression.txt"))
traj_sig_peaks_exp <- traj_sig_peaks_exp$x #675

convert_to_GRanges <- function(peak_list) {
  split_peaks <- do.call(rbind, strsplit(peak_list, "-"))
  GRanges(seqnames = split_peaks[, 1], 
          ranges = IRanges(start = as.numeric(split_peaks[, 2]), 
                           end = as.numeric(split_peaks[, 3])))
}

gr_traj <- convert_to_GRanges(traj_sig_peaks_exp)
gr_sp <- convert_to_GRanges(sp_sig_peaks)
overlaps <- findOverlaps(gr_traj, gr_sp)

overlapping_peaks_1 <- gr_traj[queryHits(overlaps)]
overlapping_peaks_2 <- gr_sp[subjectHits(overlaps)]
overlap_peaks <- data.frame(Overlap_Peaks_1 = as.character(overlapping_peaks_1), Overlap_Peaks_2 = as.character(overlapping_peaks_2)) #32
write.table(overlap_peaks, paste0(output_dir, "overlap_peaks_sp_traj_exp.txt")) 

# 3. overlap with overlapped trajectory-wise significant peaks
traj_sig_peaks_overlap <- read.table(paste0(output_dir, "significant_peaks_trajectory_overlap.txt"))
traj_sig_peaks_overlap <- gsub(":", "-", traj_sig_peaks_overlap$Overlap_Peaks_2) 

overlap_sig_peak <- intersect(sp_sig_peaks, traj_sig_peaks_overlap)
write.csv(overlap_sig_peak, paste0(output_dir, "overlap_peaks_sp_traj_overlap.csv"))





##### 9. Functional annotation of spatially significant peaks for verification
library(GenomicRanges)
library(rGREAT)
moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks.txt"))
moran_peak <- moran_result$x #1645
MERINGUE_result <- read.table(paste0(output_dir, "spatially_significant_peaks_MERINGUE.txt"))
MERINGUE_peak <- rownames(MERINGUE_result) #158
SpaGene_result <- read.table(paste0(output_dir, "spatially_significant_peaks_SpaGene.txt"))
SpaGene_peak <- rownames(SpaGene_result) #38

# spatial peaks overlap with trajectory peaks
sp_traj_overlap_peak1 <- read.table("~/spATAC_output/overlap_peaks_sp_traj_corr.txt")
sp_traj_overlap_peak1 <- sp_traj_overlap_peak1$x

# convert to gr
sp_sig_peaks <- sp_traj_overlap_peak1
gr_significant <- GRanges(
  seqnames = gsub("-.*", "", sp_sig_peaks), 
  ranges = IRanges(start = as.numeric(gsub(".*-(\\d+)-\\d+", "\\1", sp_sig_peaks)),
                   end = as.numeric(gsub(".*-\\d+-(\\d+)", "\\1", sp_sig_peaks)))
)

# run GREAT
res_BP <- great(gr_significant, "GO:BP", "TxDb.Mmusculus.UCSC.mm10.knownGene")
tb_BP <- getEnrichmentTable(res_BP)
BP_result <- tb_BP %>% select(description, genome_fraction, observed_region_hits, fold_enrichment, p_value, p_adjust) %>% filter(fold_enrichment >= 1 & p_adjust < 0.05)

# filter for pathway of interest
BP_result <- read.csv("/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/E13/spATAC_output/sp_traj_overlap_peak1_GO_result.csv")
sig_pathway <- c("growth", "central nervous system development", "developmental growth", "brain development", "negative regulation of cell differentiation", 
                 "regulation of growth", "forebrain development", "regulation of developmental growth", "cell fate commitment", "neural precursor cell proliferation", 
                 "glial cell differentiation", "cerebral cortex development", "stem cell population maintenance", "cell proliferation in forebrain", 
                 "negative regulation of neuron differentiation", "cell fate specification", "neuron fate specification", "brain development", "glial cell development",
                 "positive regulation of developmental growth", "radial glial cell differentiation", "glial cell fate specification", "negative regulation of glial cell differentiation", 
                 "cell morphogenesis involved in neuron differentiation", "neuron projection development", "positive regulation of growth", "forebrain radial glial cell differentiation", 
                 "neuron fate commitment", "forebrain cell migration")
sig_pathway_result <- BP_result %>% filter(description %in% sig_pathway)
nrow(sig_pathway_result)

# visualization
library(ggplot2)
library(forcats)
p <- ggplot(sig_pathway_result, aes(x = observed_region_hits, y = fct_reorder(description, observed_region_hits), color = fold_enrichment, size = -log10(p_adjust))) +
  geom_point(alpha = 0.6) + 
  theme_bw() +
  labs(x = "Observed Region Hits", y = "Pathway") +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", axis.text.y = element_text(size = 11.5)) +
  guides(color = guide_colorbar(title = "Fold enrichment"), size = guide_legend(title = "-log10(pval_adj)"))
ggsave("/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/E13/spATAC_output/sp_traj_overlap_peak1_GO_result.pdf", p, width = 7.5)

# BP_result <- read.csv(paste0(output_dir, "moran_peak_GO_result.csv"))
p <- ggplot(head(BP_result, 20), aes(x = observed_region_hits, y = fct_reorder(description, observed_region_hits), color = fold_enrichment, size = observed_region_hits)) +
  geom_point(alpha = 0.6) + 
  theme_bw() +
  labs(x = "Observed Region Hits", y = "Pathway", size = 12) +
  scale_color_gradient(low = "#FF9999", high = "#CC0000") +
  theme(legend.position = "right", legend.title = element_blank(), axis.text.y = element_text(size = 12)) 
ggsave(paste0(output_dir, "moran_peak_GO_result.pdf"), p)
