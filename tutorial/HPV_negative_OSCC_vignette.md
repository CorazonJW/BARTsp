# 1. Load relevant packages

HPV negative OSCC Visium raw data and processed Seurat object can be found at https://doi.org/10.6084/m9.figshare.20304456.v1. 
 
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
3. spatial_coordinates: A data frame containing the spatial locations of each cell. It must include columns named "x" and "y". 
4. cell_types: A vector of regions (cell types) of interest. 

```{r, echo=TRUE, results='markup'}
subset_object <- readRDS("~/ocss.RDS")

expression_matrix <- object@assays$SCT@counts

cell_metadata <- object@meta.data
cell_metadata$cell_type <- as.character(cell_metadata$cluster_annotations)

spatial_coordinates <- GetTissueCoordinates(object)
colnames(spatial_coordinates) <- c("x", "y")

cell_type <- c("core", "edge", "transitory")
obj <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_type)
```

# 3. Compute trajectory and obtain DEGs along pseudo-time

This step identifies genes whose expression changes as pseudo-time increases. 

```{r, echo=TRUE, results='markup'}
cds <- construct_trajectory(obj, "core")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, 
                              cor_cutoff_pos = 0.15, cor_cutoff_neg = -0.1)
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

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.1)
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
bart_proj <- bart(name = "tumor progression", genome = "mm10", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

# Increasingly-expressed genes (to predict TFs at downstream)
bart_proj <- bart(name = "tumor progression", genome = "mm10", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")
```

# 7. Visualize BART predicted results

Users can highlight the top TFs predicted by BART and/or TFs of interest in visualization step. 

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
TF_of_interest <- c("GRHL3", "TCF4", 
                    "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                    "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                    "FOS", "JUN", "JUND", "TP53")
                    
# TFs at upstream
plot_BART_results(results_geneset_up, TF_of_interest, cutoff = 0.05, top_n = 6)

# TFs at downstream
plot_BART_results(results_geneset_down, TF_of_interest, cutoff = 0.05, top_n = 6)
```

# 8. Integrate BART predicted results

This step allows integrative analysis of upstream-active and downstream-active TFs. 

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
dt <- integrate_bart_result(results_geneset_up, results_geneset_down, cutoff_up = 0.05, cutoff_down = 0.05)

# Visualize integrated ranks of all TFs
# This function shows the top n active at upstream and downstream
plot_integration_bar(dt, top_n = 10)
# This function demonstrates all TFs and highlights user-defined TFs of interst
plot_integration_dot <- (dt, tf_highlight = TF_of_interest)
```