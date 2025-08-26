
library(Seurat)

source("~/manuscript_figures.r")

subset <- readRDS("/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/E13/spRNA_output/spRNA_subset.rds")


# Figure S2D-G
DefaultAssay(subset) <- "Spatial"

p9 <- visium_plot_gene(subset, "Creb1", 50)
p10 <- visium_plot_gene(subset, "Lef1", 50)
p11 <- visium_plot_gene(subset, "Mycn", 50)
p12 <- visium_plot_gene(subset, "E2f3", 50)
ggsave("Creb1_expr.png", p9, width = 3.5, height = 4, dpi = 600)
ggsave("Lef1_expr.png", p10, width = 3.5, height = 4, dpi = 600)
ggsave("Mycn_expr.png", p11, width = 3.5, height = 4, dpi = 600)
ggsave("E2f3_expr.png", p12, width = 3.5, height = 4, dpi = 600)

# Figure S2H-K
DefaultAssay(subset) <- "peaks"

p5 <- motif_acc(subset, "CREB1")
ggsave("CREB1_motif_acc.png", p5, width = 4, height = 5, dpi = 600)
p6 <- motif_acc(subset, "LEF1")
ggsave("LEF1_motif_acc.png", p6, width = 4, height = 5, dpi = 600)
p7 <- motif_acc(subset, "MYCN")
ggsave("MYCN_motif_acc.png", p7, width = 4, height = 5, dpi = 600)
p8 <- motif_acc(subset, "E2F3")
ggsave("E2F3_motif_acc.png", p8, width = 4, height = 5, dpi = 600)


# Figure S2L-O
library(GenomicRanges)
ls <- fisher_test("CREB1", 11)
p1 <- fisher_plot(ls, "= 3.477e-10", 2)
ggsave("CREB1_motif_fisher.png", p1, width = 2.8, height = 2, dpi = 600)
ls <- fisher_test("LEF1", 14)
p2 <- fisher_plot(ls, "= 6.335e-05", 4)
ggsave("LEF1_motif_fisher.png", p2, width = 2.8, height = 2, dpi = 600)
ls <- fisher_test("MYCN", 8)
p3 <- fisher_plot(ls, "< 2.2e-16", 6)
ggsave("MYCN_motif_fisher.png", p3, width = 2.8, height = 2, dpi = 600)
ls <- fisher_test("E2F3", 10)
p4 <- fisher_plot(ls, "< 2.2e-16", 20)
ggsave("E2F3_motif_fisher.png", p4, width = 2.8, height = 2, dpi = 600)
