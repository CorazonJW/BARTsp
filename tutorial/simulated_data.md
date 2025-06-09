# 1. Load relevant packages

Simulated data can be found at https://github.com/CorazonJW/BARTsp/tree/main/data/simulated_data/.
 
```{r, echo=TRUE, results='markup'}
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

library(BARTsp)
```

# 2. Prepare input

```{r, echo=TRUE, results='markup'}
expression_matrix <- read.table("expression_matrix.txt")
expression_matrix <- as.matrix(expression_matrix)
expression_matrix <- as(expression_matrix, "dgCMatrix")

cell_metadata <- read.table("cell_metadata.txt")

feature_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

spatial_coordinates <- read.table("spatial_coordinates.txt")

obj <- prepare_input(expression_matrix, cell_metadata, feature_metadata, spatial_coordinates, cell_types = unique(cell_metadata$cell_type))
```

# 3. Compute trajectory and obtain DEGs along pseudo-time

```{r, echo=TRUE, results='markup'}
cds <- construct_trajectory(obj, "1")
pseudotime_values <- monocle3::pseudotime(cds)

traj_DEG <- get_traj_features(pseudotime_values, obj, pval_cutoff = 0.1, 
                              cor_cutoff_pos = 0.25, cor_cutoff_neg = -0.15)
```

# 4. Compute Moran's I for each gene and select DEGs across spatial locations

```{r, echo=TRUE, results='markup'}
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

moran_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.15)
```

# 5. Construct input for BART

```{r, echo=TRUE, results='markup'}
genes <- get_sig_features_geneset(traj_DEG, moran_DEG)
geneset <- construct_BART_geneset_input(genes)
```

# 6. Run BART

Decreasingly-expressed genes are used to predict TFs active at upstream in the trajectory
Increasingly-expressed genes are used to predict TFs active at downstream in the trajectory

```{r, echo=TRUE, results='markup'}
# Decreasingly-expressed genes (to predict TFs at upstream )
bart_proj <- bart(name = "GRN", genome = "hg38", data = geneset$down_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_up <- get_BART_results(bart_proj, "geneset")

# Increasingly-expressed genes (to predict TFs at downstream)
bart_proj <- bart(name = "GRN", genome = "hg38", data = geneset$up_gene$significant_features, type = "geneset")
bart_proj <- run_BART(bart_proj, type = "geneset")
results_geneset_down <- get_BART_results(bart_proj, "geneset")
```

# 7. Visualize BART predicted results

```{r, echo=TRUE, results='markup', fig.width=10, fig.height=8}
# TFs at upstream
plot_BART_results(results_geneset_up, TF_of_interest = NULL, cutoff = 0.05, top_n = 6)

# TFs at downstream
plot_BART_results(results_geneset_down, TF_of_interest = NULL, cutoff = 0.05, top_n = 6)
```
