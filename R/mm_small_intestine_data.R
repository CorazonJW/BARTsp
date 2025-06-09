#' Mouse Small Intestine Visium HD Dataset
#'
#' A subset of the mouse small intestine Visium HD dataset containing cell type labels.
#' This dataset is derived from the 10x Genomics Visium HD CytAssist gene expression
#' libraries of mouse intestine.
#'
#' @format A Seurat object with the following components:
#' \describe{
#'   \item{assays}{Contains the spatial gene expression data}
#'   \item{meta.data}{Contains cell metadata including cell type labels}
#'   \item{coordinates}{Contains spatial coordinates for each spot}
#' }
#'
#' @source \url{https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine}
#'
#' @examples
#' data(mm_small_intestine_data)
#' # View cell type distribution
#' table(mm_small_intestine_data$enterocyte_type)
"mm_small_intestine_data" 