# BARTsp

BARTsp is an R package for analyzing spatial transcriptomics/epigenomics data with BART integration. It provides tools for identifying spatially variable genes, constructing trajectories, and predicting transcription factors.

## Installation

You can install the development version of BARTsp from GitHub:

```r
# install.packages("devtools")
devtools::install_github("CorazonJW/BARTsp")
```

## Features

- Spatial gene expression analysis using SPARKX, Moran's I, or KNN
- Trajectory inference using Monocle3
- BART integration for transcription factor prediction
- Visualization tools for results

## Usage

Here's a basic example of how to use BARTsp:

```r
library(BARTsp)

# Input (Visium V2)
expression_matrix <- seu_object@assays$Spatial.008um@layers$counts
colnames(expression_matrix) <- rownames(seu_object@assays$Spatial.008um@cells)
rownames(expression_matrix) <- rownames(seu_object@assays$Spatial.008um@features)
cell_metadata <- seu_object@meta.data # cell_metadata must contain a column named "Cell_type"
spatial_coordinates <- GetTissueCoordinates(seu_object)

# Prepare input data
input_data <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_types)

# Run trajectory analysis
trajectory <- construct_trajectory(input_data, start_cell_type = "cell_type_A")
traj_DEG <- get_traj_features(trajectory$pseudotime, input_data, pval_cutoff = 0.05, 
                              cor_cutoff_pos = 0.3, cor_cutoff_neg = -0.3)

# Run spatial analysis (Moran's I)
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

sp_DEG <- get_moran_result(morana_I_result, adj.val = 0.1, moransI = 0.15)

# (Alternative) Run spatial analysis (SPARKX)
sparkx_result <- run_SPARKX(input_data, numCores = 4)
sp_DEG <- get_sparkx_DEGs(sparkx_result, cutoff = 0.05)

# (Alternative) Run spatial analysis (KNN)
knn_result <- run_knn_spatial(input_data, k = 5, method = "correlation", cutoff = 0.5)

# Construct BART input
genes <- get_sig_features_geneset(traj_DEG, sp_DEG)
bart_input <- construct_BART_geneset_input(genes)

# Run BART analysis
bart_proj <- bart(name = "my_analysis", genome = "mm10", data = bart_input, type = "geneset")
bart_results <- run_BART(bart_proj)
results <- get_BART_results(bart_proj)

# Visualize results
plot_BART_results(results, TF_of_interest = c("TF1", "TF2"), cutoff = 0.1)

# Integrate BART predicted results
dt <- integrate_bart_result(results_geneset_up, results_geneset_down, cutoff_up = 0.05, cutoff_down = 0.05)

plot_integration_bar(dt, top_n = 10)
plot_integration_dot <- (dt, tf_highlight = TF_of_interest)
```

## Documentation

For detailed documentation, please visit the [package website](https://CorazonJW.github.io/BARTsp/).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 