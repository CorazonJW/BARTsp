
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ArchR)
library(grid)

setwd("~/spATAC_output/")
data_dir <- "~/data/"
output_dir <- "~/spATAC_output/"

##### 1. Pre-process
# Load data
threads = 8
addArchRThreads(threads = threads)
addArchRGenome("mm10")
inputFiles <- paste0(data_dir, "peak_expression_matrix.tsv.gz")
sampleNames <- 'ME13'
# Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000), 
  force = TRUE
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE, 
  force = TRUE
)
# Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data_dir, "ME13_50um_spatial"), filter.matrix = filter.matrix)
name <- paste0("ME13#", rownames(image@coordinates), "-1")
rownames(image@coordinates) <- name
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue
# Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)
proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)
proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
proj_in_tissue <- addImputeWeights(proj_in_tissue)



##### 2. Add MOCA predicted id to ArchR object from processed Seurat object
E13_sp <- readRDS("~/spRNA_output/spRNA_preprocess.rds")

dt <- data.frame(E13_sp$predicted.id)
rownames(dt) <- paste0("ME13#", rownames(dt), "-1")
dt <- dt[match(proj_in_tissue@cellColData@rownames, rownames(dt)), , drop = FALSE]

proj_in_tissue$predicted.id <- dt$E13_sp.predicted.id



##### 3. Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "predicted.id")

pathToMacs2 <- findMacs2()

proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "predicted.id", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

saveRDS(proj_in_tissue, paste0(output_dir, "spatac_preprocess.RDS"))
proj_in_tissue <- readRDS(paste0(output_dir, "spatac_preprocess.RDS"))



###### 4. Motif analysis to add motif information to ArchR object
# ChromVAR Deviatons Enrichment
if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
    proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif")
}

proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveRDS(proj_in_tissue, paste0(output_dir, "spatac_preprocess.RDS"))
proj_in_tissue <- readRDS(paste0(output_dir, "spatac_preprocess.RDS"))



##### 5. Subset project to only contain radial glia and postmitotic premature neurons
E13_spRNA <- readRDS("~/spRNA_output/spRNA_preprocess.rds")
Idents(E13_spRNA) <- E13_spRNA$predicted.id
subset <- subset(E13_spRNA, idents = c("Radial glia", "Postmitotic premature neurons"))
cell_subset <- paste0("ME13#", rownames(subset@meta.data), "-1")
ct_label <- as.character(subset$predicted.id)

subset_ATAC <- subsetArchRProject(
  ArchRProj = proj_in_tissue,
  cells = cell_subset,
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

saveRDS(subset_ATAC, paste0(output_dir, "spatac_subset.RDS"))
subset_ATAC <- readRDS(paste0(output_dir, "spatac_subset.RDS"))





##### 6. Construct trajectory
subset_ATAC <- readRDS(paste0(output_dir, "spatac_subset.RDS"))
trajectory <- c("Radial glia", "Postmitotic premature neurons")

proj_traj <- addTrajectory(
  ArchRProj = subset_ATAC, 
  name = "trajectory", 
  groupBy = "predicted.id",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
head(proj_traj$trajectory[!is.na(proj_traj$trajectory)])

saveRDS(proj_traj, paste0(output_dir, "spatac_trajectory.RDS"))
proj_traj <- readRDS(paste0(output_dir, "spatac_trajectory.RDS"))

# Plot Motif
plotVarDev <- getVarDeviations(proj_traj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_traj, addDOC = FALSE)

# Visualize spATAC trajectory 
trajMM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajPM  <- getTrajectory(ArchRProj = proj_traj, name = "trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)

# spatac_gene <- sub("chr[0-9XY]+:", "", trajGSM@NAMES)
# spatac_motif <- gsub("^(deviations:|z:)|_[0-9]+$", "", trajMM@NAMES)
# spatac_peak <- gsub("chr(\\d+):(\\d+)_(\\d+)", "chr1-\\2-\\3", trajPM@NAMES)

pdf(paste0(output_dir, "pseudotime_motif.pdf"), width = 10, height = 8)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
p1
p2
p3
dev.off()

# Visualize spATAC trajectory 
source("~/manuscript_figures.r")

subset <- readRDS("~/spRNA_output/spRNA_subset.rds")

dt <- data.frame(cell = proj_traj@cellColData@rownames, sp_ATAC_traj = proj_traj$trajectory)
dt$cell <- gsub("^ME13#|-1$", "", dt$cell)
dt <- dt[match(WhichCells(subset), dt$cell), ]

subset <- AddMetaData(object = subset, metadata = dt$sp_ATAC_traj, col.name = "spATAC_traj")
DefaultAssay(subset) <- "ATAC"
subset <- subset[, !is.na(subset$spATAC_traj)]

saveRDS(subset, "~/spRNA_subset.rds")

p4 <- visium_ptime_spatial(subset, 18, "spATAC_traj")
ggsave("pseudotime.png", p4, width = 5.5, height = 4, dpi = 600)

