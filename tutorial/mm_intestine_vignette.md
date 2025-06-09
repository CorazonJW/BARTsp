# 1. Load relevant packages

Mouse small intestine Visium HD raw data can be found at https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine. Processed Seurat object can be found at https://github.com/CorazonJW/BARTsp/tree/main/data/mm_small_intestine_data.RDS.
 
```{r, echo=TRUE, results='markup'}
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

library(BARTsp)
```

# 2. Prepare input

Inputs include
1. expression_matrix: A gene by cell read count matrix. 
2. cell_metadata: A data frame of cell information. It must contain a column named "cell_type". 
3. feature_metadata: A data frame of feature names, required by Monocle3 for pseudo-time analysis. 
4. spatial_coordinates: A data frame containing the spatial locations of each cell. It must include columns named "x" and "y". 
5. cell_types: A vector of regions (cell types) of interest. 

```{r, echo=TRUE, results='markup'}
subset_object <- readRDS("~/mm_small_intestine.RDS")

expression_matrix <- object@assays$Spatial.008um@layers$counts
colnames(expression_matrix) <- rownames(object@assays$Spatial.008um@cells)
rownames(expression_matrix) <- rownames(object@assays$Spatial.008um@features)
cell_metadata <- object@meta.data
cell_metadata$cell_type <- cell_metadata$enterocyte_type
feature_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
spatial_coordinates <- GetTissueCoordinates(object)

obj <- prepare_input(expression_matrix, cell_metadata, feature_metadata, spatial_coordinates, cell_types = c("Enterocyte_Progenitor", "Enterocyte_Immature", "Enterocyte_Mature"))
```

# 3. Compute trajectory and obtain DEGs along pseudo-time

This step identifies genes whose expression changes as pseudo-time increases. 

```{r, echo=TRUE, results='markup'}
cds <- construct_trajectory(obj, start_cell_type = "Enterocyte_Progenitor")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, 
                              cor_cutoff_pos = 0.05, cor_cutoff_neg = -0.05)
```

# 4. Compute Moran's I for each gene and select DEGs across spatial locations

This step calculates Moran's I for each gene. Moran's I is a measure of spatial autocorrelation. It tells you whether similar values (e.g., gene expression) tend to cluster together in space. A positive Moran's I indicates that the gene is spatially variable, meaning its expression is not random but clustered in certain regions. This spatial variability reflects underlying cell differentiation patterns, where specific genes are upregulated in localized populations of differentiating cells.

```{r, echo=TRUE, results='markup'}
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.01)
```

## (Alternative methods to find SVGs)

Users can use SPARKX or KNN-based methods to identify spatially variable genes. 

```{r, echo=TRUE, results='markup'}
# SPARKX
sparkx_result <- run_SPARKX(obj, numCores = 4)
sparkx_DEG <- get_sparkx_DEGs(sparkx_result, cutoff = 0.05)

# KNN-based method
knn_result <- run_knn_spatial(obj, k = 5, method = "correlation", cutoff = 0.3)
```

# 5. Construct input for BART

In geneset mode, we use the intersection of pseudo-time-related DEGs and SVGs to ensure precision. The overlapping gene set is then divided into two categories based on their expression patterns along the trajectory. Genes whose expression increases with pseudo-time are considered downstream-active genes, likely regulated by transcription factors active later in the differentiation process. In contrast, genes whose expression negatively correlates with pseudo-time are considered upstream-active genes, assumed to be regulated by TFs acting earlier in the trajectory.

```{r, echo=TRUE, results='markup'}
genes <- get_sig_features_geneset(traj_DEG, moran_DEG)
geneset <- construct_BART_geneset_input(genes)
```

# 6. Run BART

Decreasingly-expressed genes are used to predict TFs active at upstream in the trajectory
Increasingly-expressed genes are used to predict TFs active at downstream in the trajectory

```{r, echo=TRUE, results='markup'}
# Decreasingly-expressed genes (to predict TFs at upstream )
bart_proj <- bart(name = "enterocyte", genome = "mm10", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

# Increasingly-expressed genes (to predict TFs at downstream)
bart_proj <- bart(name = "radial glia to PPN", genome = "mm10", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")
```

# 7. Visualize BART predicted results

Users can highlight the top TFs predicted by BART and/or TFs of interest in visualization step. 

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
TF_of_interest <- c("HNF4G", "HNF4A", "GATA6", "HNF1B", "MAF", "CDX2", "MAFB", "GATA4", "ATOH1", "HES1")

# TFs at upstream
plot_BART_results(results_geneset_up, TF_of_interest, cutoff = 0.05, top_n = 6)

# TFs at downstream
plot_BART_results(results_geneset_down, TF_of_interest, cutoff = 0.05, top_n = 6)
```
