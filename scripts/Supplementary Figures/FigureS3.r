
##### Codes for generating Figure S3D-L can be found under ~/scripts/Figure6/Compare_tools/for_SVFs_detection.r

spATAC_bart_result_rep1 <- read.csv("~/rep1/BART_results_upstream.csv") %>% mutate(rank = X)
spATAC_bart_result_rep2 <- read.csv("~/rep2/BART_results_upstream.csv") %>% mutate(rank = X)

tf_highlight <-  c("PAX6", "SOX9", "FEZF2", "KLF4", "NEUROD2")
p1 <- plot_BART_results(spATAC_bart_result_rep1, tf_highlight, 0.05, 6)
p2 <- plot_BART_results(spATAC_bart_result_rep2, tf_highlight, 0.05, 6)

ggsave("BART_ATAC_rep1_prediction.png", p1, width = 3, height = 3, dpi = 600)
ggsave("BART_ATAC_rep2_prediction.png", p2, width = 3, height = 3, dpi = 600)


subset_seurat_obj <- readRDS("~/spatac_subset_seurat_obj_rep1.RDS")

Idents(subset_seurat_obj) <- subset_seurat_obj$cell_type

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

p1 <- visium_cell_type_tissue(subset_seurat_obj, pallete, 5, "bottom")
ggsave("tissue_spatial_map.png", p1, width = 8, height = 4, dpi = 600)
