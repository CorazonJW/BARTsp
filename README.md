# BART-spatial

BART-spatial is an R package for detecting functional transcription regulators based on spatial transcriptomics/epigenomics data. It integrates spatial variability and pseudo-temporal information with publicly available TF binding profiles, boosting prediction accuracy in the absence of high TF expression. 
Applied to multiple datasets, BART-spatial outperformed existing methods, identifying stage-specific TFs and revealing regulators undetectable by expression alone. Its compatibility with spatial epigenomics data further strengthens its prediction power and enables cross-validation. 

## Installation

You can install BART-spatial from GitHub:

```r
# install.packages("devtools")
devtools::install_github("CorazonJW/BARTsp")
```

## Features

- Spatial gene expression analysis using Moran's I
- Pseudo-time inference using Monocle3
- BART integration for transcription factor prediction

## Usage

Here's a basic example of how to use BARTsp:

### 1. Prepare inputs
```r
library(BARTsp)

# Three inputs are required: (1) expression matrix; (2) cell metadata; (3) spatial coordinates
seurat_object <- readRDS("PATH/TO/SEURAT_OBJECT")
expression_matrix <- seurat_object@assays$Spatial.008um@layers$counts
colnames(expression_matrix) <- rownames(seurat_object@assays$Spatial.008um@cells)
rownames(expression_matrix) <- rownames(seurat_object@assays$Spatial.008um@features)
cell_metadata <- seurat_object@meta.data # cell_metadata must contain a column named "Cell_type"
spatial_coordinates <- GetTissueCoordinates(seurat_object)

input_data <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_types)
```

### 2. Detect pseudo-temporally varaible features (TVFs)
```r
trajectory <- construct_trajectory(input_data, start_cell_type = "cell_type_A")
traj_DEG <- get_traj_features(trajectory$pseudotime, input_data, pval_cutoff = 0.05, 
                              cor_cutoff_pos = 0.1, cor_cutoff_neg = -0.1)
```

### 3. Detect spatially varaible features (SVFs)
```r
# Moran's I
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

sp_DEG <- get_moran_result(morana_I_result, adj.val = 0.05, moransI = 0.1)

# (Alternative) SPARKX
sparkx_result <- run_SPARKX(input_data, numCores = 4)
sp_DEG <- get_sparkx_DEGs(sparkx_result, cutoff = 0.05)

# (Alternative) KNN
knn_result <- run_knn_spatial(input_data, k = 5, method = "correlation", cutoff = 0.1)
```

### 4. Construct BART algorithm input
```r
genes <- get_sig_features_geneset(traj_DEG, sp_DEG)
bart_input <- construct_BART_geneset_input(genes)
```

### 5. Run BART
```r
bart_proj <- bart(name = "my_analysis", genome = "mm10", data = bart_input, type = "geneset")
bart_results <- run_BART(bart_proj)
results <- get_BART_results(bart_proj)

plot_BART_results(results, TF_of_interest = c("TF1", "TF2"), cutoff = 0.05)
```

### 6. Integrate upstream and downstream predictions
```r
dt <- integrate_bart_result(results_geneset_up, results_geneset_down, cutoff_up = 0.05, cutoff_down = 0.05)

plot_integration_bar(dt, top_n = 10)
plot_integration_dot <- (dt, tf_highlight = TF_of_interest)
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 