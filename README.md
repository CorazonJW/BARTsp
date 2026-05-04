# BART-spatial

BART-spatial is an R package for detecting functional transcription regulators based on spatial transcriptomics/epigenomics data. It integrates spatial variability and pseudo-temporal information with publicly available TF binding profiles, boosting prediction accuracy in the absence of high TF expression. 
Applied to multiple datasets, BART-spatial outperformed existing methods, identifying stage-specific TFs and revealing regulators undetectable by expression alone. Its compatibility with spatial epigenomics data further strengthens its prediction power and enables cross-validation. 

## Installation

Before installing BART-spatial, please install BARTsc. More details about BARTsc can be seen under https://github.com/hongpan-uva/BARTsc. 

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

## Usage

Please refer to [tutorial](./tutorial) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 