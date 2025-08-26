
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
## Load the normalized data 
E13_sp <- readRDS(paste0(data_dir, "/E13_spatial_RNA_ATAC.rds"))
E13_sp <- FindVariableFeatures(E13_sp, selection.method = "vst")

## QC
p1 <- SpatialFeaturePlot(E13_sp, features = "nCount_Spatial", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + 
            guides(fill = guide_colorbar(title = "Read Count")) +
            theme(legend.title = element_text(size = 14, face = "bold"), 
                  legend.text = element_text(size = 10), 
                  legend.key.width = unit(1, "cm"), 
                  legend.key.height = unit(0.5, "cm"))
p2 <- SpatialFeaturePlot(E13_sp, features = "nFeature_Spatial", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + 
            guides(fill = guide_colorbar(title = "Feature Count")) +
            theme(legend.title = element_text(size = 14, face = "bold"), 
                  legend.text = element_text(size = 10), 
                  legend.key.width = unit(1, "cm"), 
                  legend.key.height = unit(0.5, "cm"))
p <- p1+p2
ggsave(paste0(output_dir, "nCount_nFeature.pdf"), p, width = 8, height = 5)

Idents(E13_sp) <- E13_sp$orig.ident
p1 <- VlnPlot(E13_sp, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(E13_sp, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p <- p1+p2
ggsave(paste0(output_dir, "nCount_nFeature_violin_plot.pdf"), p, width = 4, height = 8)

E13_sp <- subset(E13_sp, nCount_Spatial > 200 & nCount_Spatial < 15000 & nFeature_Spatial > 200 & nFeature_Spatial < 5000)



##### 2. Define cell types by MOCA projection
## Prepare MOCA reference
MOCA_count <- readRDS(paste0(data_dir, "MOCA_dt/gene_count_cleaned.RDS"))
MOCA_gene_name <- read.table(paste0(data_dir, "MOCA_dt/gene_annotate.csv"), header = T, sep = ",") %>% select(gene_id, gene_short_name)
rownames(MOCA_count) <- MOCA_gene_name$gene_short_name[match(rownames(MOCA_count), MOCA_gene_name$gene_id)]
MOCA_ct_label <- read.table(paste0(data_dir, "MOCA_dt/cell_annotate.csv"), header = T, sep = ",") %>% select(sample, Main_cell_type, development_stage)
MOCA_ct_label <- MOCA_ct_label[MOCA_ct_label$development_stage == 13.5, ]
# clean data
MOCA_count <- MOCA_count[!duplicated(rownames(MOCA_count)), ] # keep the first occurance
rownames(MOCA_count) <- gsub("_", "-", rownames(MOCA_count))
rownames(MOCA_count) <- trimws(rownames(MOCA_count)) # remove leading/trailing whitespace
MOCA_count <- MOCA_count[, colnames(MOCA_count) %in% MOCA_ct_label$sample]
# create Seurat object
MOCA_obj <- CreateSeuratObject(counts = MOCA_count)
MOCA_obj$cell_type <- MOCA_ct_label$Main_cell_type[match(colnames(MOCA_obj), MOCA_ct_label$sample)]
MOCA_obj <- NormalizeData(MOCA_obj)
MOCA_obj <- FindVariableFeatures(MOCA_obj, selection.method = "vst")

## QC
Idents(MOCA_obj) <- MOCA_obj$orig.ident
p1 <- VlnPlot(MOCA_obj, features = "nCount_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p2 <- VlnPlot(MOCA_obj, features = "nFeature_RNA", pt.size = 0.1, alpha = 0) + NoLegend()
p <- p1+p2
ggsave("~/data/nCount_nFeature_violin_plot.pdf", p, width = 4, height = 8)

MOCA_obj_filter <- subset(MOCA_obj, nFeature_RNA > 300 & nFeature_RNA < 1500)

saveRDS(MOCA_obj_filter, "~/data/MOCA_dt/E13_sp_MOCA_data_filtered.rds")


## Projection
options(future.globals.maxSize = 10 * 1024^3) 
MOCA_obj <- readRDS(paste0(data_dir, "MOCA_dt/E13_sp_MOCA_data_filtered.rds"))
MOCA_obj@meta.data <- MOCA_obj@meta.data[!is.na(MOCA_obj@meta.data$cell_type), ]

DefaultAssay(E13_sp) <- "Spatial"
anchors <- FindTransferAnchors(reference = MOCA_obj, query = E13_sp) 
transferred_label <- TransferData(anchorset = anchors, refdata = MOCA_obj$cell_type) 
E13_sp <- AddMetaData(object = E13_sp, metadata = transferred_label$predicted.id, col.name = "predicted.id")
E13_sp <- AddMetaData(object = E13_sp, metadata = transferred_label$prediction.score.max, col.name = "predicted.score")

saveRDS(E13_sp, paste0(output_dir, "spRNA_preprocess.rds"))


## Visualization
source("~/manuscript_figures.r")

E13_sp <- readRDS("~/spRNA_preprocess.rds")
Idents(E13_sp) <- E13_sp$predicted.id

object <- subset(E13_sp, idents = c("Radial glia", "Postmitotic premature neurons"), invert = FALSE)
Idents(object) <- object$predicted.id

# Prediction score
p <- SpatialFeaturePlot(E13_sp, features = "predicted.score", alpha = c(0.9, 0.9), pt.size.factor = 1.4) + 
            guides(fill = guide_colorbar(title = "Prediction score")) +
            theme(legend.title = element_text(size = 14, face = "bold"), 
                  legend.text = element_text(size = 10), 
                  legend.key.width = unit(1, "cm"), 
                  legend.key.height = unit(0.5, "cm"))
ggsave(paste0(output_dir, "spatial_ct_prediction_score.pdf"), p, width = 8, height = 5)

Idents(E13_sp) <- E13_sp$predicted.id
p <- VlnPlot(E13_sp, features = "predicted.score", pt.size = 0.1) + NoLegend()
ggsave(paste0(output_dir, "ct_prediction_score.pdf"), p, width = 8, height = 5)

# Cell type labels
pallete <- c("Stromal cells"="#E6194B", "Osteoblasts"="#3CB44B", "Myocytes"="#FFE119", 
             "Connective tissue progenitors"="#4363D8", "Excitatory neurons"="#F58231", 
             "Radial glia"="#E5CCFF", "Postmitotic premature neurons"="#4D1A89", 
             "Definitive erythroid lineage"="#46F0F0", "Chondroctye progenitors"="#BCF60C", 
             "Endothelial cells"="#FABEBE", "Notochord cells"="#008080", "Neural progenitor cells"="#E6BEFF", 
             "Cholinergic neurons"="#9A6324", "Epithelial cells"="#FFFAC8", "Oligodendrocyte Progenitors"="#800000", 
             "Isthmic organizer cells"="#AA6E28", "Jaw and tooth progenitors"="#808000", "Sensory neurons"="#FFD8B1", 
             "Inhibitory neurons"="#000075", "Granule neurons"="#FF007F", "Ependymal cell"="#FCCDE5", 
             "Inhibitory interneurons"="#9933FF", "Premature oligodendrocyte"="#1B9E77", 
             "Inhibitory neuron progenitors"="#D95F02")

p1 <- visium_cell_type_tissue(E13_sp, pallete, 45, "right")
ggsave("tissue_spatial_map.png", p1, width = 8, height = 4, dpi = 600)

p2 <-  visium_cell_type_spatial(object, pallete, 45)
ggsave("spatial_map.png", p2, width = 5.5, height = 4, dpi = 600)

# Marker gene expression
p4 <- visium_plot_gene(object, "Pax6", 50)
p5 <- visium_plot_gene(object, "Fabp7", 50)
p6 <- visium_plot_gene(object, "Myt1l", 50)
ggsave("Pax6_expr.png", p4, width = 3.5, height = 4, dpi = 600)
ggsave("Fabp7_expr.png", p5, width = 3.5, height = 4, dpi = 600)
ggsave("Myt1l_expr.png", p6, width = 3.5, height = 4, dpi = 600)
