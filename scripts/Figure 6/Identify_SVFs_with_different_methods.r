

##### Method 0: Test with SPARK (not working)
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
# got error, cannot resolve





##### Method 1: Test with Seurat Moran's I
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





##### Method 2: Test with Moran's I to find spatially differential expressed peaks
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
significant_result <- morana_I_result[sapply(morana_I_result, function(x) x$adjusted_p.value < 0.05)] #581

sp_sig_peaks <- names(significant_result)
write.table(sp_sig_peaks, paste0(output_dir, "spatially_significant_peaks.txt"))

positive_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] > x$estimate[[2]])] # all positive regions
negative_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] < x$estimate[[2]])]





##### Method 3: Test with singleCellHaystack to find spatially differential expressed peaks (not very well)
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





##### Method 4: Test with SPADE to find spatially differential expressed peaks (not very well)
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





##### 2. Check overlap of MERINGUE results, Seurat Moran's I test, and self-written Moran's I
seurat_moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks_seurat_moran's_I.txt"))
seurat_moran_peak <- rownames(seurat_moran_result) #2166
moran_result <- read.table(paste0(output_dir, "spatially_significant_peaks.txt"))
moran_peak <- moran_result$x #581
MERINGUE_result <- read.table(paste0(output_dir, "spatially_significant_peaks_MERINGUE.txt"))
MERINGUE_peak <- rownames(MERINGUE_result) #158

overlap_peak <- intersect(seurat_moran_peak, moran_peak) #507 (good overlap)
overlap_peak <- intersect(MERINGUE_peak, seurat_moran_peak) #140 (good overlap)
