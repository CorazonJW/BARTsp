
# packages required for Visium HD
# install.packages("hdf5r")
# install.packages("arrow")
# remotes::install_version(package = 'Seurat', version = package_version('5.1.0'))

library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(Matrix)

data_dir <- "~/data/"
output_dir <- "~/results/"

options(future.globals.maxSize = 20 * 1024^3)  # 15 GB

##### 1. Pre-process
## Normalization
object <- Load10X_Spatial(data.dir = data_dir, bin.size = c(2, 8, 16))

Assays(object)

for(Res in c("002","008","016")){
    ASname <- paste0("Spatial.", Res, "um")
    DefaultAssay(object) <- ASname
    object <- NormalizeData(object)
}

saveRDS(object, file=paste0(output_dir, "small_intestine_preprocess.rds"))

sampleName <- rep(c("002um","008um","016um"))
sampleCol <- c("#FF5555","#D43F00","#8B0000")

pdf(file=paste0(output_dir, "sampleQC.pdf"), width=6, height=6)
par(mfrow=c(1,3),mar=c(8,4,2,2))
barplot(log10(c(length(which(!is.na(object$nCount_Spatial.002um))), 
                 length(which(!is.na(object$nCount_Spatial.008um))), 
                 length(which(!is.na(object$nCount_Spatial.016um))))),
        names.arg = sampleName, ylab = "log10(#spot)", main = "#valid spot", col = sampleCol, las = 2)
boxplot(log10(as.numeric(object$nCount_Spatial.002um)+1), 
        log10(as.numeric(object$nCount_Spatial.008um)+1),
        log10(as.numeric(object$nCount_Spatial.016um)+1),
        names=sampleName,ylab="log10 reads count",main="#reads per spot",col=sampleCol, outline=F,las=2)
boxplot(log10(as.numeric(object$nFeature_Spatial.002um)+1), 
        log10(as.numeric(object$nFeature_Spatial.008um)+1),
        log10(as.numeric(object$nFeature_Spatial.016um)+1),
        names=sampleName,ylab="log10 covered gene",main="#genes per spot",col=sampleCol, outline=F,las=2)
abline(h=log10(c(10,100,500)), lwd=2)
dev.off()


## QC
# nFeature >= 10 for 002um; nFeature >= 100 for 008um
object <- readRDS(paste0(output_dir, "small_intestine_preprocess.rds"))

DefaultAssay(object) <- "Spatial.002um"
object_002um_highQ <- object[, which(object$nFeature_Spatial.002um >=10)]
DefaultAssay(object) <- "Spatial.008um"
object_008um_highQ <- object[, which(object$nFeature_Spatial.008um >=100)]
DefaultAssay(object) <- "Spatial.016um"
object_016um_highQ <- object[, which(object$nFeature_Spatial.016um >=500)]

## Preprocess
preprocess <- function(inobject, assay){
  object <- FindVariableFeatures(inobject)
  object <- ScaleData(object)
  object <- RunPCA(object, assay = as.character(assay), reduction.name = "pca.spatial")
  object <- FindNeighbors(object, assay = as.character(assay), reduction = "pca.spatial", dims = 1:50)
  object <- RunUMAP(object, reduction = "pca.spatial", reduction.name = "umap.spatial", return.model = T, dims = 1:50)
  for(RES in seq(0.5, 3, 0.5)){
    object <- FindClusters(object, resolution = RES)
  }
  return(object)
}
# object_002um_highQ <- preprocess(object_002um_highQ, assay = "Spatial.002um")
object_008um_highQ <- preprocess(object_008um_highQ, assay = "Spatial.008um")
object_016um_highQ <- preprocess(object_016um_highQ, assay = "Spatial.016um")

saveRDS(object_002um_highQ, file=paste0(output_dir, "object_002um_highQ.rds"))
saveRDS(object_008um_highQ, file=paste0(output_dir,"object_008um_highQ.rds"))
saveRDS(object_016um_highQ, file=paste0(output_dir,"object_016um_highQ.rds"))


object_002um_highQ <- readRDS(paste0(output_dir,"object_002um_highQ.rds"))
object_008um_highQ <- readRDS(paste0(output_dir,"object_008um_highQ.rds"))
object_016um_highQ <- readRDS(paste0(output_dir,"object_016um_highQ.rds"))

# q02 <-  FeaturePlot(object_002um_highQ, label = F,features="nFeature_Spatial.008um", reduction = "umap.spatial",cols = c("lightgrey", "darkgreen")) + ggtitle("002um #covGene")
q08 <- FeaturePlot(object_008um_highQ, label = F,features="nFeature_Spatial.008um", reduction = "umap.spatial",cols = c("lightgrey", "darkgreen")) + ggtitle("008um #covGene")
q16 <- FeaturePlot(object_016um_highQ, label = F,features="nFeature_Spatial.016um", reduction = "umap.spatial",cols = c("lightgrey", "darkgreen")) + ggtitle("016um #covGene")
png(file=paste0(output_dir, "UMAP_features_covGene_p1.png"),width=1000,height=500)
grid.arrange(q08,q16,ncol=2,nrow=1)
dev.off()


p08_1 <- DimPlot(object_008um_highQ, label = TRUE,group.by="Spatial.008um_snn_res.0.5", reduction = "umap.spatial") + ggtitle("008um res0.5")
p08_2 <- DimPlot(object_008um_highQ, label = TRUE,group.by="Spatial.008um_snn_res.1", reduction = "umap.spatial") + ggtitle("008um res1.0")
p08_3 <- DimPlot(object_008um_highQ, label = TRUE,group.by="Spatial.008um_snn_res.1.5", reduction = "umap.spatial") + ggtitle("008um res1.5")
p08_4 <- DimPlot(object_008um_highQ, label = TRUE,group.by="Spatial.008um_snn_res.2", reduction = "umap.spatial") + ggtitle("008um res2.0")

png(file=paste0(output_dir, "UMAP_features_clusterRes_008um.png"),width=1000,height=1000)
grid.arrange(p08_1,p08_2,p08_3,p08_4,ncol=2,nrow=2)
dev.off()

p16_1 <- DimPlot(object_016um_highQ, label = TRUE,group.by="Spatial.016um_snn_res.0.5", reduction = "umap.spatial") + ggtitle("016um res0.5")
p16_2 <- DimPlot(object_016um_highQ, label = TRUE,group.by="Spatial.016um_snn_res.1", reduction = "umap.spatial") + ggtitle("016um res1.0")
p16_3 <- DimPlot(object_016um_highQ, label = TRUE,group.by="Spatial.016um_snn_res.1.5", reduction = "umap.spatial") + ggtitle("016um res1.5")
p16_4 <- DimPlot(object_016um_highQ, label = TRUE,group.by="Spatial.016um_snn_res.2", reduction = "umap.spatial") + ggtitle("016um res2.0")

png(file=paste0(output_dir, "UMAP_features_clusterRes_016um.png"),width=1000,height=1000)
grid.arrange(p16_1,p16_2,p16_3,p16_4,ncol=2,nrow=2)
dev.off()


Idents(object_008um_highQ) <- "Spatial.008um_snn_res.0.5"
DefaultAssay(object_008um_highQ) <- "Spatial.008um"
p1 <- DimPlot(object_008um_highQ, reduction = "umap.spatial", label = F) + theme(legend.position = "bottom")
p2 <- SpatialDimPlot(object_008um_highQ, label = F) + theme(legend.position = "bottom")
p <- p1 | p2
ggsave(paste0(output_dir, "clusters_on_image_008um.png"), p)

Idents(object_016um_highQ) <- "Spatial.016um_snn_res.0.5"
DefaultAssay(object_016um_highQ) <- "Spatial.016um"
p1 <- DimPlot(object_016um_highQ, reduction = "umap.spatial", label = F) + theme(legend.position = "bottom")
p2 <- SpatialDimPlot(object_016um_highQ, label = F) + theme(legend.position = "bottom")
p <- p1 | p2
ggsave(paste0(output_dir, "clusters_on_image_016um.png"), p)




##### 2. Cell type label annotation
object_008um_highQ <- readRDS("~/results/object_008um_highQ.rds")

## Create Seurat object of reference scRNA-seq data
data_path <- "~/data/sc_ref/GSE92332_atlas_UMIcounts.txt.gz"
umi_matrix <- read.delim(gzfile(data_path), row.names = 1, check.names = FALSE)
umi_sparse <- as(as.matrix(umi_matrix), "dgCMatrix")

seurat_obj <- CreateSeuratObject(counts = umi_sparse, project = "GSE92332", min.cells = 3, min.features = 200)

cell_types <- sapply(strsplit(colnames(seurat_obj), "_"), function(x) tail(x, 1))
seurat_obj$cell_type <- cell_types
head(seurat_obj@meta.data$cell_type)
table(seurat_obj$cell_type)

saveRDS(seurat_obj, "~/data/sc_ref/mm_intestine_ref.RDS")

ref <- readRDS("~/data/sc_ref/mm_intestine_ref.RDS")
ref <- SCTransform(ref, verbose = TRUE)
ref <- RunPCA(ref, assay = "SCT", verbose = TRUE)

## Projection
ref.anchors <- FindTransferAnchors(reference = ref, query = object_008um_highQ, dims = 1:20, reference.reduction = "pca", normalization.method = "SCT")
predicted_ct <- TransferData(anchorset = ref.anchors, refdata = ref$cell_type, dims = 1:20)

object_008um_highQ <- AddMetaData(object_008um_highQ, metadata = predicted_ct$predicted.id, col.name = "Cell_type")
object_008um_highQ <- AddMetaData(object_008um_highQ, metadata = predicted_ct$prediction.score.max, col.name = "Prediction_score")

saveRDS(object_008um_highQ, paste0(output_dir, "object_008um_highQ_with_ct_label.RDS"))


## Visualization
p <- ggplot(object_008um_highQ@meta.data, aes(x = Prediction_score)) +
        geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
        theme_minimal()   
ggsave(paste0(output_dir, "prediction_score.pdf"), p)

Idents(object_008um_highQ) <- object_008um_highQ$Cell_type
p <- SpatialDimPlot(object_008um_highQ, label = F) + theme(legend.position = "bottom")
ggsave(paste0(output_dir, "cell_type.pdf"), p, width = 10)


coords <- GetTissueCoordinates(object_008um_highQ)
p <- ggplot(coords, aes(x = x, y = y)) + geom_point(size = 0.1) + coord_fixed() + theme_bw()
ggsave(paste0(output_dir, "coords.pdf"), p, width = 8, height = 6)

selected_barcodes <- rownames(coords)[coords$y >= 10000 & coords$y <= 12500 & coords$x >= 13500 & coords$x <= 15000]

subset_object <- subset(object_008um_highQ, cells = selected_barcodes)
pallette <- c("Enterocyte.Immature.Distal"  = "#E41A1C", "Goblet" = "#377EB8", "Enterocyte.Mature.Distal" = "#4DAF4A", 
              "Stem" = "#984EA3", "Paneth" = "#FF7F00", "Enterocyte.Progenitor.Early" = "#FFFF33", 
              "Enterocyte.Progenitor.Late" = "#A65628", "Endocrine" = "#F781BF", "TA.Early" = "#009E73", 
              "Enterocyte.Progenitor" = "#66C2A5", "TA.G1" = "#FC8D62", "Tuft" = "#8DA0CB", "TA.G2" = "#FFD92F")

Idents(subset_object) <- subset_object$Cell_type
p <- SpatialDimPlot(subset_object, label = F, pt.size.factor = 18, image.alpha = 0.7, alpha = c(0.8, 0.8), cols = pallette) + 
        labs(fill = "Cell Type") +
        theme(legend.text = element_text(size = 12))
ggsave(paste0(output_dir, "cell_type_susbet.pdf"), p, width = 8, height = 6)

saveRDS(subset_object, paste0(output_dir, "subset_object_with_ct_label.RDS"))



subset_ct_object <- subset(subset_object, subset = Cell_type %in% c("Enterocyte.Progenitor", "Enterocyte.Progenitor.Late", "Enterocyte.Progenitor.Early", 
                                                                    "Enterocyte.Immature.Distal", "Enterocyte.Mature.Distal"))

subset_ct_object$enterocyte_type <- ifelse(subset_ct_object$Cell_type %in% c("Enterocyte.Progenitor", "Enterocyte.Progenitor.Late", "Enterocyte.Progenitor.Early"), "Enterocyte_Progenitor",
                                        ifelse(subset_ct_object$Cell_type == "Enterocyte.Immature.Distal", "Enterocyte_Immature", 
                                            ifelse(subset_ct_object$Cell_type == "Enterocyte.Mature.Distal", "Enterocyte_Mature", NA)))

pallette <- c("Enterocyte_Immature"  = "#E41A1C", "Enterocyte_Mature" = "#4DAF4A", "Enterocyte_Progenitor" = "#FFFF33")

Idents(subset_ct_object) <- subset_ct_object$enterocyte_type
p <- SpatialDimPlot(subset_ct_object, label = F, pt.size.factor = 18, image.alpha = 0.7, alpha = c(0.8, 0.8), cols = pallette) + 
        labs(fill = "Cell Type") +
        theme(legend.text = element_text(size = 12))
ggsave(paste0(output_dir, "subset_enterocyte.pdf"), p, width = 8, height = 6)
