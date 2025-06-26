
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)

setwd("~/spRNA_output/")
data_dir <- "~/data/"
output_dir <- "~/spRNA_output/"

##### 1. Pre-process
# Load the processed data 
E13_sp <- readRDS(paste0(data_dir, "/E13_spatial_RNA_ATAC.rds"))
E13_sp <- FindVariableFeatures(E13_sp, selection.method = "vst")

# E13_sp <- subset(E13_sp, nCount_Spatial > 200 & nCount_Spatial < 15000 & nFeature_Spatial > 200 & nFeature_Spatial < 5000)

# visualize nFeature and nCount
p1 <- SpatialFeaturePlot(E13_sp, features = "nCount_Spatial", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + guides(fill = guide_colorbar(title = "Read Count")) +
            theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.5, "cm"))
p2 <- SpatialFeaturePlot(E13_sp, features = "nFeature_Spatial", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + guides(fill = guide_colorbar(title = "Feature Count")) +
            theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.5, "cm"))
p <- p1+p2
ggsave(paste0(output_dir, "nCount_nFeature.pdf"), p, width = 8, height = 5)

Idents(E13_sp) <- E13_sp$orig.ident
p1 <- VlnPlot(E13_sp, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(E13_sp, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p <- p1+p2
ggsave(paste0(output_dir, "nCount_nFeature_violin_plot.pdf"), p, width = 4, height = 8)




##### 2. Define cell types by MOCA projection
# Prepare MOCA reference
MOCA_count <- readRDS(paste0(data_dir, "MOCA_dt/gene_count_cleaned.RDS"))
# change gene id to name
MOCA_gene_name <- read.table(paste0(data_dir, "MOCA_dt/gene_annotate.csv"), header = T, sep = ",") %>% select(gene_id, gene_short_name)
rownames(MOCA_count) <- MOCA_gene_name$gene_short_name[match(rownames(MOCA_count), MOCA_gene_name$gene_id)]
# import cell type annotation
MOCA_ct_label <- read.table(paste0(data_dir, "MOCA_dt/cell_annotate.csv"), header = T, sep = ",") %>% select(sample, Main_cell_type, development_stage)
MOCA_ct_label <- MOCA_ct_label[MOCA_ct_label$development_stage == 13.5, ]
# remove duplicates
MOCA_count <- MOCA_count[!duplicated(rownames(MOCA_count)), ] # keep the first occurance
rownames(MOCA_count) <- gsub("_", "-", rownames(MOCA_count)) # replace underscores with dashes
rownames(MOCA_count) <- trimws(rownames(MOCA_count)) # remove leading/trailing whitespace
MOCA_count <- MOCA_count[, colnames(MOCA_count) %in% MOCA_ct_label$sample]
# create Seurat object
MOCA_obj <- CreateSeuratObject(counts = MOCA_count)
# add cell type label
MOCA_obj$cell_type <- MOCA_ct_label$Main_cell_type[match(colnames(MOCA_obj), MOCA_ct_label$sample)]
# process MOCA reference
MOCA_obj <- NormalizeData(MOCA_obj)
MOCA_obj <- FindVariableFeatures(MOCA_obj, selection.method = "vst")

Idents(MOCA_obj) <- MOCA_obj$orig.ident
p1 <- VlnPlot(MOCA_obj, features = "nCount_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p2 <- VlnPlot(MOCA_obj, features = "nFeature_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p <- p1+p2
ggsave("~/data/nCount_nFeature_violin_plot.pdf", p, width = 4, height = 8)

MOCA_obj_filter <- subset(MOCA_obj, nFeature_RNA > 300 & nFeature_RNA < 1500)

# visualize nFeature and nCount
Idents(MOCA_obj_filter) <- MOCA_obj_filter$orig.ident
p1 <- VlnPlot(MOCA_obj_filter, features = "nCount_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p2 <- VlnPlot(MOCA_obj_filter, features = "nFeature_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p <- p1+p2
ggsave("~/data/MOCA_dt/nCount_nFeature_violin_plot_after_QC.pdf", p, width = 4, height = 8)

saveRDS(MOCA_obj_filter, "~/data/MOCA_dt/E13_sp_MOCA_data_filtered.rds")


# Projection
options(future.globals.maxSize = 10 * 1024^3)  # 10GB
MOCA_obj <- readRDS(paste0(data_dir, "MOCA_dt/E13_sp_MOCA_data_filtered.rds"))
MOCA_obj@meta.data <- MOCA_obj@meta.data[!is.na(MOCA_obj@meta.data$cell_type), ]

DefaultAssay(E13_sp) <- "Spatial"
anchors <- FindTransferAnchors(reference = MOCA_obj, query = E13_sp) # find anchor
transferred_label <- TransferData(anchorset = anchors, refdata = MOCA_obj$cell_type) # transfer vague labels
# transferred_label <- transferred_label %>% mutate(new_predict.id = ifelse(prediction.score.max < 0.5, "unknown", predicted.id))
E13_sp <- AddMetaData(object = E13_sp, metadata = transferred_label$predicted.id, col.name = "predicted.id")
E13_sp <- AddMetaData(object = E13_sp, metadata = transferred_label$prediction.score.max, col.name = "predicted.score")

saveRDS(E13_sp, paste0(output_dir, "spRNA_preprocess.rds"))
E13_sp <- readRDS(paste0(output_dir, "spRNA_preprocess.rds"))

# visualize prediction score 
p <- SpatialFeaturePlot(E13_sp, features = "predicted.score", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + guides(fill = guide_colorbar(title = "Prediction score")) +
            theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.5, "cm"))
ggsave(paste0(output_dir, "spatial_ct_prediction_score.pdf"), p, width = 8, height = 5)

Idents(E13_sp) <- E13_sp$predicted.id
p <- VlnPlot(E13_sp, features = "predicted.score", pt.size = 0.1) + NoLegend()
ggsave(paste0(output_dir, "ct_prediction_score.pdf"), p, width = 8, height = 5)

# visualize predicted cell type
pallete <- c(
  "Stromal cells"="#E6194B", "Osteoblasts"="#3CB44B", "Myocytes"="#FFE119", "Connective tissue progenitors"="#4363D8", "Excitatory neurons"="#F58231", "Radial glia"="#911EB4", 
  "Definitive erythroid lineage"="#46F0F0", "Postmitotic premature neurons"="#F032E6", "Chondroctye progenitors"="#BCF60C", "Endothelial cells"="#FABEBE", "Notochord cells"="#008080", "Neural progenitor cells"="#E6BEFF", 
  "Cholinergic neurons"="#9A6324", "Epithelial cells"="#FFFAC8", "Oligodendrocyte Progenitors"="#800000", "Isthmic organizer cells"="#AA6E28", "Jaw and tooth progenitors"="#808000", "Sensory neurons"="#FFD8B1", 
  "Inhibitory neurons"="#000075", "Granule neurons"="#FF007F", "Ependymal cell"="#FCCDE5", "Inhibitory interneurons"="#9933FF", "Premature oligodendrocyte"="#1B9E77", "Inhibitory neuron progenitors"="#D95F02"
)

p1 <- SpatialDimPlot(E13_sp, group.by = "Jiont_clusters", alpha = c(0.8, 0.8), pt.size.factor = 1.5) + labs(fill = "Joint cluster")
p2 <- SpatialDimPlot(E13_sp, group.by = "predicted.id", cols = pallete, alpha = c(0.8, 0.8), pt.size.factor = 1.5) + labs(fill = "Cell type")
p <- p1+p2
ggsave(paste0(output_dir, "cell-type_label.pdf"), p, width = 16, height = 8)

# visualize QC metrics of each cell type
Idents(E13_sp) <- E13_sp$predicted.id
p1 <- VlnPlot(E13_sp, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(E13_sp, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p <- p1+p2
ggsave(paste0(output_dir, "nCount_nFeature_violin_plot_ct.pdf"), p, width = 8, height = 10)






##### 3. Trajectory analysis to find genes related to radial glia --> Postmitotic prematrue neurons differentiation
# Subset raidal glia and Postmitotic prematrue neurons (radial glia --> Postmitotic prematrue neurons)
Idents(E13_sp) <- E13_sp$predicted.id
subset <- subset(E13_sp, idents = c("Radial glia", "Postmitotic premature neurons"), invert = FALSE)

p1 <- SpatialDimPlot(subset, group.by = "Jiont_clusters", alpha = c(0.8, 0.8), pt.size.factor = 1.5) + labs(fill = "Joint cluster")
p2 <- SpatialDimPlot(subset, group.by = "predicted.id", cols = pallete, alpha = c(0.8, 0.8), pt.size.factor = 1.5) + labs(fill = "Cell type")
p <- p1+p2
ggsave(paste0(output_dir, "cell-type_label_subset.pdf"), p, width = 16, height = 8)

# Check known marker gene expression across spatial locations
for (i in c("Pax6", "Fabp7", "Sox9", "Myt1l", "Tubb3", "Bcan", "Luzp2")) {
    p <- SpatialFeaturePlot(subset, features = i, pt.size.factor = 1.5)
            theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.2, "cm"))
    ggsave(paste0(output_dir, "DEG_expression/", i, "_exp.pdf"), p, width = 5, height = 4)
}

saveRDS(subset, paste0(output_dir, "spRNA_subset.rds"))
subset <- readRDS(paste0(output_dir, "spRNA_subset.rds"))

# Derive pseudo-time by monocle3 (summarized as 'construct_trajectory' function in BARTsp)
library(SeuratWrappers)
library(monocle3)

# prepare inputs
expression_matrix <- subset@assays$Spatial@counts
cell_metadata <- subset@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
# create monocle3 object
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
# preprocess
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
# learn trajectory graph
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
# construct trajectory
get_earliest_principal_node <- function(cds, time_bin="Radial glia"){
  cell_ids <- which(colData(cds)[, "predicted.id"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
# Plot trajectory
p1 <- plot_cells(cds, color_cells_by = "cluster")
p2 <- plot_cells(cds, color_cells_by="predicted.id")
p3 <- plot_cells(cds, color_cells_by = "pseudotime")
p <- p1+p2+p3
ggsave(paste0(output_dir, "pseudotime-analysis.pdf"), p, width = 12, height = 4)
# extract pseudotime
pseudotime <- pseudotime(cds)
subset$pseudotime <- pseudotime[colnames(subset)]

# Visualize on spatial map
p <- SpatialFeaturePlot(subset, features = "pseudotime")
ggsave(paste0(output_dir, "pseudotime.pdf"), p, width = 5, height = 4)

saveRDS(subset, paste0(output_dir, "spRNA_subset.rds"))
subset <- readRDS(paste0(output_dir, "spRNA_subset.rds"))

# Find variable genes along pseudotime by correlation analysis (summarized as 'get_traj_features' function in BARTsp)
expression_matrix <- subset@assays$Spatial@counts
pseudotime_values <- subset$pseudotime
cor_results <- apply(expression_matrix, 1, function(gene_expr) {
    cor.test(gene_expr, pseudotime_values, method = "spearman")$p.value
})
adjusted_pvals <- p.adjust(cor_results, method = "fdr")
significant_genes <- names(adjusted_pvals[adjusted_pvals < 0.01]) 
significant_genes <- significant_genes[!is.na(significant_genes)]

write.table(significant_genes, paste0(output_dir, "significant_genes_trajectory.txt"))







##### 4. Spatailly variable genes analysis to find genes whose expression varies in radial glia and PPN region
# Use Moran'I to test SVGs (summarized as 'prepare_moran_input', 'preprocess_data', and 'compute_morans_I' functions in BARTsp)
library(Seurat)
library(spdep)
library(Matrix)
library(expm)
library(stats)
library(dplyr)
library(parallel)

# Load data ('prepare_moran_input')
# expression matrix
subset <- readRDS("~/spRNA_output/spRNA_subset.rds")
gene <- subset@assays$Spatial@counts
genes <- rownames(gene) 
# spatial coordinates
spatial_loc <- read.table("~/data/Spots.coords.tsv") %>% filter(V1 %in% colnames(subset))
colnames(spatial_loc) <- c("cell", "x", "y")
rownames(spatial_loc) <- spatial_loc$cell
spatial_loc <- as.data.frame(spatial_loc)
spatial_loc$x <- as.numeric(spatial_loc$x)
spatial_loc$y <- as.numeric(spatial_loc$y)
S <- as.matrix(spatial_loc[, -1])

# Center and scale gene expression and spatial coordinate matrix ('preprocess_data')
preprocess_data <- function(expression_matrix, spatial_coordinates) {
  y <- scale(expression_matrix, center = TRUE, scale = TRUE)
  S <- scale(spatial_coordinates, center = TRUE, scale = TRUE)

  return(list(y = y, S = S))
}

# Compute Moran's I for spatial autocorrelation 'compute_morans_I'
compute_morans_I <- function(y, S) {
  coords <- as.data.frame(S)
  k <- min(5, nrow(coords) - 1)
  neighbors <- knearneigh(coords, k = k)  # Find nearest neighbors
  listw <- nb2listw(knn2nb(neighbors), style = "W")
  
  morans_test <- moran.test(y, listw)
  return(morans_test)
}

# Process all genes 
morana_I_result <- list()
for (i in genes) {
    y <- as.matrix(gene[as.character(i), ])
    colnames(y) <- as.character(i)
    # Skip if all values in y are zero
    if (all(y == 0, na.rm = TRUE)) {
        message(paste("Skipping gene:", i, " - all values are zero"))
        next
    }
    # Skip if y contains NA values
    if (any(is.na(y))) {
        message(paste("Skipping gene:", i, " - contains NA values"))
        next
    }
    processed_data <- preprocess_data(y, S)
    # Skip if NA values remain after processing
    if (any(is.na(processed_data$y))) {
        message(paste("Skipping gene:", i, "- NA values in scaled expression"))
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

# Select significant genes < (summarized as 'get_moran_result' function)
significant_result <- morana_I_result[sapply(morana_I_result, function(x) x$estimate[[1]] != x$estimate[[2]])] #15772/16048
significant_result <- morana_I_result[sapply(morana_I_result, function(x) x$adjusted_p.value < 0.05)] #3702
sp_sig_genes <- names(significant_result)
write.table(sp_sig_genes, paste0(output_dir, "SVGs_moran's_I.txt"))

positive_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] > x$estimate[[2]])] # all positive regions
negative_result <- significant_result[sapply(significant_result, function(x) x$estimate[[1]] < x$estimate[[2]])]