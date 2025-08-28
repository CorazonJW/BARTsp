# BART-spatial

BART-spatial is an R package for detecting functional transcription regulators based on spatial transcriptomics/epigenomics data. It integrates spatial variability and pseudo-temporal information with publicly available TF binding profiles, boosting prediction accuracy in the absence of high TF expression. 
Applied to multiple datasets, BART-spatial outperformed existing methods, identifying stage-specific TFs and revealing regulators undetectable by expression alone. Its compatibility with spatial epigenomics data further strengthens its prediction power and enables cross-validation. 

## Installation

You can install BART-spatial from GitHub:

```r
# Install BARTsc first
install.packages("devtools")

devtools::install_github("immunogenomics/presto")
devtools::install_github("hongpan-uva/BARTsc")
```

After BARTsc is installed, for the first time the user imports it, BARTsc needs to be initialized with function initialize(). This step will automatically create a python virtual environment and install the BART2 python module and related dependencies. The user can specify the path for storing relevant data library (recommended) and the path for the module.

```r
library("BARTsc")

initialize()
```

```r
# Then install BART-spatial
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

# Three inputs are required: (1) expression matrix; (2) cell metadata containing a column named "Cell_type"; (3) spatial coordinates

input_data <- prepare_input(expression_matrix, cell_metadata, spatial_coordinates, cell_types)
```

### 2. Detect pseudo-temporally varaible features (TVFs)
```r
# For RNA mode, monocle3 is used to compute pseudo-time. 
trajectory <- construct_trajectory(input_data, start_cell_type = "cell_type_A")

# For ATAC mode, we recommand users to use ArchR to calculate pseudo-time. 

TVFs <- get_traj_features(trajectory$pseudotime, input_data, pval_cutoff = 0.05, 
                          cor_cutoff_pos = 0.1, cor_cutoff_neg = -0.1)
```

### 3. Detect spatially varaible features (SVFs)
```r
moran_obj <- prepare_moran_input(obj)
moran_obj <- preprocess_data(moran_obj)
morana_I_result <- compute_morans_I(moran_obj)

p_values <- sapply(morana_I_result, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "fdr")
for (i in seq_along(morana_I_result)) {
    morana_I_result[[i]]$adjusted_p.value <- adjusted_p_values[i]
}

SVFs <- get_moran_result(morana_I_result, adj.val = 0.05, moransI = 0.1)
```

### 4. Construct BART algorithm input
#### RNA mode
```r
genes <- get_sig_features_geneset(TVFs, SVFs)
bart_input <- construct_BART_geneset_input(genes)
```
#### ATAC mode
```r
region <- get_sig_features_region(obj, TVFs, SVFs)
bart_input <- construct_BART_region_input(region, description_up = "upstream", description_down = "downstream")
```

### 5. Run BART
#### RNA mode
```r
bart_proj <- bart(name = "my_analysis", genome = "mm10", data = bart_input, type = "geneset")
bart_results <- run_BART(bart_proj, type = "geneset")
results <- get_BART_results(bart_proj, type = "geneset")

plot_BART_results(results, TF_of_interest = c("TF1", "TF2"), cutoff = 0.05)
```
#### ATAC mode
```r
bart_proj <- bart(name = "my_analysis", genome = "mm10", data = bart_input, type = "region")
bart_results <- run_BART(bart_proj, type = "region")
results <- get_BART_results(bart_proj, type = "region")

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