## 1. Load relevant packages

E13 mouse embryo co-profiling data can be found at https://cells.ucsc.edu/?ds=brain-spatial-omics+e13. Single-cell RNA reference can be found at https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads. 
 
```{r, echo=TRUE, results='markup'}
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

library(BARTsp)
```

## 2. Prepare input

Inputs include
1. expression_matrix: A gene by cell read count matrix. 
2. cell_metadata: A data frame of cell information. It must contain a column named "cell_type". 
3. spatial_coordinates: A data frame containing the spatial locations of each cell. It must include columns named "x" and "y". 
4. cell_types: A vector of regions (cell types) of interest. 

```{r, echo=TRUE, results='markup'}
E13_sp <- readRDS("~/mm_embryo_spRNA_spATAC.RDS")

expression_matrix <- E13_sp@assays$Spatial@counts

cell_metadata <- E13_sp@meta.data
cell_metadata$cell_type <- cell_metadata$predicted.id

spatial_coordinates <- data.frame(E13_sp@images$slice1@coordinates)
spatial_coordinates$x <- spatial_coordinates$imagerow
spatial_coordinates$y <- spatial_coordinates$imagecol

obj <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, 
                     c("Radial glia", "Postmitotic premature neurons"))
```

## 3. Detect pseudo-temporally variable features (TVFs)

This step identifies genes whose expression changes as pseudo-time increases. 

```{r, echo=TRUE, results='markup'}
cds <- construct_trajectory(obj, "Radial glia")
pseudotime_values <- monocle3::pseudotime(cds)

TVFs <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, 
                          cor_cutoff_pos = 0.1250701, cor_cutoff_neg = -0.2456667)
```

## 4. Detect spatially variable features (SVFs)

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

SVFs <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.2373)
```

## 5. Construct input for BART algorithm

In geneset mode, we use the intersection of TVFs and SVFs to ensure precision. The overlapping geneset is then divided into two categories based on their expression patterns along the trajectory. Genes whose expression increases with pseudo-time are considered downstream-active genes, likely regulated by transcription regulators (TRs) active later in the differentiation process. In contrast, genes whose expression negatively correlates with pseudo-time are considered upstream-active genes, assumed to be regulated by TRs acting earlier in the trajectory.

```{r, echo=TRUE, results='markup'}
genes <- get_sig_features_geneset(TVFs, SVFs)
geneset <- construct_BART_geneset_input(genes)
```

## 6. Run BART

Decreasingly-expressed genes are used to predict TRs active at upstream in the trajectory. 
Increasingly-expressed genes are used to predict TRs active at downstream in the trajectory. 

```{r, echo=TRUE, results='markup'}
# Decreasingly-expressed genes (to predict TRs at upstream)
bart_proj <- bart(name = "Upstream", genome = "mm10", data = geneset$upstream$feature, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

# Increasingly-expressed genes (to predict TRs at downstream)
bart_proj <- bart(name = "Downstream", genome = "mm10", data = geneset$downstream$feature, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")
```

## 7. Visualize BART-spatial prediction

Users can highlight the top TRs predicted by BART and/or TRs of interest in visualization step. 

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
TF_of_interest <- c("PAX6", "SOX9", "NEUROD2", "KLF4", "FEZF2","HES1")

# TRs at upstream
plot_BART_results(results_geneset_up, TF_of_interest, 0.05, 6)

# TRs at downstream
plot_BART_results(results_geneset_up, TF_of_interest, 0.05, 6)
```

## 8. Integrate BART-spatial prediction

This step allows integrative analysis of upstream-active and downstream-active TRs. 

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
dt <- integrate_bart_result(results_geneset_up, results_geneset_down, cutoff_up = 0.05, cutoff_down = 0.05)

# This function shows the top n active at upstream and downstream
plot_integration_bar(dt, top_n = 10)

# This function demonstrates all upstream and downsteam predictions as a single plot and highlights user-defined TRs of interst
plot_integration_dot(dt, tf_highlight = TF_of_interest)
```