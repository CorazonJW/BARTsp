library(BARTsp)
library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)
library(sciColors)
library(forcats)
library(stringr)
library(gridExtra)

##### Visium data visualization
cell_type_spatial <- function(cell_df, x, y, cell_type, pallette) {
    ggplot(cell_df, aes(x = {{x}}, y = {{y}}, color = {{cell_type}})) +
        geom_point(size = 2) +
        scale_color_manual(values = pallette) +
        labs(color = "Stage", x = "", y = "") +
        theme_bw() +  
        theme(legend.position = "top", 
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(0,0,0,0),
              legend.title = element_text(size = 18, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 16, face = "bold", family = "Arial"),
              legend.key.width = unit(0.5, "cm"), 
              legend.spacing.x = unit(0.1, "cm"),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank()) +
        guides(color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE))
}

ptime_spatial <- function(cell_df, x, y, pseudotime_values) {
    ggplot(cell_df, aes(x = {{x}}, y = {{y}}, color = {{pseudotime_values}})) +
        geom_point(size = 2) +
        scale_color_viridis_c(option = "magma", name = "Pseudotime") +
        labs(color = "Pseudotime", x = "", y = "") +
        theme_bw() +  
        theme(legend.position = "top", 
              legend.key.width = unit(1, "cm"), 
              legend.key.height = unit(0.5, "cm"),
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(0,0,0,0),
              legend.title = element_text(size = 18, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 16, face = "bold", family = "Arial"), 
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank()) +
        guides(color = guide_colorbar(title.position = "top", barwidth = unit(6, "cm"), 
                                      barheight = unit(0.4, "cm")))
}

plot_gene <- function(expression_matrix, gene_name, cell_df, x, y) {
  expr <- expression_matrix[gene_name, ]
  df <- cell_df
  df$expr <- expr
  p <- ggplot(cell_df, aes(x = {{x}}, y = {{y}}, color = expr)) +
        geom_point(size = 2) +
        scale_color_gradient(low = "#E0E0E0", high = "#FF0000") +
        labs(color = paste0(gene_name,"\nexpr"), x = "", y = "") + 
        theme_bw() +
        theme(legend.position = "top", 
              legend.key.width = unit(0.8, "cm"), 
              legend.key.height = unit(0.5, "cm"),
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(0,0,0,0),
              legend.title = element_text(size = 20, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 18, face = "bold", family = "Arial"), 
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
  return(p)
}

avg_exp_input_gene <- function(cell_df, x, y, input_gene, title) {
    ggplot(cell_df, aes(x = {{x}}, y = {{y}}, color = {{input_gene}})) +
        geom_point(size = 2.5) +
        scale_color_viridis_c() +
        labs(color = title, x = "", y = "") +
        theme_bw() +
        theme(legend.position = "top", 
              legend.key.width = unit(0.85, "cm"), 
              legend.key.height = unit(0.5, "cm"),
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(0,0,0,0),
              legend.title = element_text(size = 16, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 14, family = "Arial"),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
}












##### Visium HD data visualization
visium_cell_type_tissue <- function(seu_object, pallette, size, legend) {
    SpatialDimPlot(seu_object, label = F, pt.size.factor = size, image.alpha = 0.8, alpha = c(0.8, 0.8), cols = pallette) + 
        labs(fill = "Cell Type") + 
        theme(legend.position = legend, 
              legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 12, family = "Arial"),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
}

visium_cell_type_spatial <- function(seu_object, pallette, size) {
    SpatialDimPlot(seu_object, label = F, pt.size.factor = as.numeric(size), image.alpha = 0.8, alpha = c(0.8, 0.8), cols = pallette) + 
        labs(fill = "Cell Type") + 
        theme(legend.position = "right", 
              legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 12, family = "Arial"), 
              legend.key.height = unit(0.5, "cm"), 
              legend.spacing.y = unit(0.15, "cm"),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
}

visium_ptime_spatial <- function(seu_object, size, ptime, labs) {
    SpatialFeaturePlot(seu_object, as.character(ptime), pt.size.factor = as.numeric(size), image.alpha = 0.8, alpha = c(0.8, 0.8)) +
        theme(legend.position = "top", 
              legend.key.width = unit(1, "cm"), 
              legend.key.height = unit(0.5, "cm"),
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(0,0,0,0),
              legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 12, family = "Arial"), 
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank()) +
        labs(fill = labs)
}

visium_plot_gene <- function(seu_object, gene_name, size) {
    expr_values <- FetchData(seu_object, vars = gene_name, slot = "data")[[1]]
    upper_lim <- quantile(expr_values, 0.99, na.rm = TRUE)

    SpatialFeaturePlot(seu_object, features = gene_name, pt.size.factor = as.numeric(size), slot = "data") + 
        labs(fill = paste0(gene_name, "  \nexpr  "), x = "", y = "") + 
        guides(fill = guide_colorbar(direction = "horizontal", title.position = "left",
               barwidth = unit(5, "cm"), barheight = unit(0.5, "cm"), ticks = FALSE)) +
        scale_fill_gradientn(colors = c("#E0E0E0", "#FF0000"), 
                             limits = c(0, upper_lim), oob = scales::squish, guide = "colorbar") +
        theme(legend.position = "top", 
              legend.key.width = unit(2, "cm"), 
              legend.key.height = unit(0.8, "cm"),
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(2,2,2,2),
              legend.title = element_text(size = 18, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 16, family = "Arial"), 
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
    # ggsave(paste0(output_dir, gene_name, "_exp.png"), p, width = 5, height = 4)
}














##### Benchmarking plots
VS_barplot <- function(dt){
    ggplot(dt, aes(x = method, y = n_TFs, fill = method)) +
            geom_col(width = 0.7) +
            geom_text(aes(label = in_TF), vjust = 1.3, size = 5, color = "#202020") +
            scale_fill_manual(values = palette, drop = FALSE) + 
            labs(x = "", y = "Number of functional TFs") +
            theme_bw() +
            theme(legend.position = "none", 
                  axis.text = element_text(size = 14, family = "Arial"), 
                  axis.title.x = element_blank(), 
                  axis.title.y = element_text(size = 18, face = "bold", family = "Arial"))
}


scale_rank <- function(dt, tf_interest) {
    dt$TF_type <- ifelse(dt$TF %in% tf_interest, "Functional", "Other prediction")
    dt$TF_type <- factor(dt$TF_type, levels = c("Functional", "Other prediction"))
    dt$rank <- as.numeric(dt$rank)
    dt$scaled_rank <- round((max(dt$rank) - dt$rank) / (max(dt$rank) - min(dt$rank)), 3)

    return(dt)
}


VS_boxplot <- function(combined_rank) {
    combined_rank$TF_type <- factor(combined_rank$TF_type, levels = c("Functional", "Other prediction"))

    plot <- ggplot(combined_rank, aes(x = method, y = scaled_rank, fill = TF_type)) +
            geom_boxplot(position = position_dodge(width = 0.7), width = 0.6) +
            scale_fill_manual(values = c("Functional" = "#E64B35", "Other prediction" = "gray70")) +
            labs(x = "", y = "Scaled rank", fill = "TF class") +
            theme_bw() +
            theme(legend.position = "right", 
                  legend.title = element_text(size = 16, face = "bold", family = "Arial"), 
                  legend.text = element_text(size = 14, family = "Arial"), 
                  axis.text = element_text(size = 14, family = "Arial"), 
                  axis.title.x = element_blank(), 
                  axis.title.y = element_text(size = 18, face = "bold", family = "Arial"))

    annotation_label <- function(p) {
    if (p < 0.00001) "*****"
    else if (p < 0.0001) "****"
    else if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else "NS"
    }

    stats <- combined_rank %>%
                group_by(method) %>%
                summarize(p = wilcox.test(scaled_rank ~ TF_type, exact = FALSE, alternative = "greater")$p.value,
                          y = max(scaled_rank, na.rm = TRUE) - 0.01, .groups = "drop") %>%
                mutate(star = vapply(p, annotation_label, character(1)))

    p_stat <- plot + geom_text(data = stats, aes(x = method, y = y, label = star), vjust = 0, size = 3, inherit.aes = FALSE)

    return(p_stat)
}

















##### Fisher's exact test of TF motifs
fisher_test <- function(TF_name, length_overlap){
    # UDHS
    mm10_UDHS <- read.table("/project/zanglab_project/jw4xtu/BART_library/mm10_library/bart2_mm10_UDHS.bed")
    mm10_UDHS_gr <- GRanges(seqnames = mm10_UDHS$V1, ranges = IRanges(start = mm10_UDHS$V2 + 1, end = mm10_UDHS$V3))

    # spatially variable peaks on UDHS
    moran_DAR <- read.csv("/project/zanglab_project/jw4xtu/spatial_project/spatial-ATAC-RNA-seq/final_results/spATAC/sp_regionset.csv")

    sp_peaks <- moran_DAR$significant_features
    peak_info <- do.call(rbind, strsplit(sp_peaks, "-"))
    output_df <- data.frame(V1 = peak_info[,1], V2 = as.integer(peak_info[,2]), V3 = as.integer(peak_info[,3]))
    rownames(output_df) <- paste0(output_df$V1, "-", output_df$V2, "-", output_df$V3)
    sp_peaks_gr <- GRanges(seqnames = output_df$V1, ranges = IRanges(start = output_df$V2, end = output_df$V3))
    sp_peaks_on_UDHS <- findOverlaps(sp_peaks_gr, mm10_UDHS_gr, minoverlap = 50)
    sp_overlap_indices <- subjectHits(sp_peaks_on_UDHS)
    sp_peaks_UDHS <- mm10_UDHS_gr[sp_overlap_indices]

    # spatially non-variable peaks on UDHS
    non_sp_overlap_indices <- setdiff(seq_along(mm10_UDHS_gr), sp_overlap_indices)
    non_sp_peaks_UDHS <- mm10_UDHS_gr[non_sp_overlap_indices]

    bed_path <- paste0("/project/zanglab_project/jw4xtu/motif_sites/mappable_results/", TF_name, "_map.bed")
    bed_files <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(bed_files) <- c("chr", "start", "end", "seq", "pval", "strand")
    motif_sites <- GRanges(seqnames = bed_files$chr, ranges = IRanges(start = bed_files$start + 1, end = bed_files$end)) # 298391
    motif_sites_on_UDHS <- findOverlaps(mm10_UDHS_gr, motif_sites, minoverlap = as.numeric(length_overlap))
    motif_sites_overlap_indices <- queryHits(motif_sites_on_UDHS)
    motif_sites_UDHS <- mm10_UDHS_gr[motif_sites_overlap_indices]

    ### Fisher's exact test
    sp_with_TF <- findOverlaps(motif_sites_UDHS, sp_peaks_UDHS, type = "within")
    nonsp_with_TF <- findOverlaps(motif_sites_UDHS, non_sp_peaks_UDHS, type = "within") 

    a <- length(sp_with_TF)
    b <- length(sp_peaks_UDHS) - length(sp_with_TF)
    c <- length(nonsp_with_TF)
    d <- length(non_sp_peaks_UDHS) - length(nonsp_with_TF)

    contingency_table <- matrix(c(a, b, c, d), nrow = 2, dimnames = list(Spatial = c("KLF4+", "KLF4-"), PeakType = c("Spatial", "Non-Spatial")))
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    print(fisher_result)

    return(list(sp_with_TF = sp_with_TF, nonsp_with_TF = nonsp_with_TF, sp_peaks_UDHS = sp_peaks_UDHS, non_sp_peaks_UDHS = non_sp_peaks_UDHS))
}

fisher_plot <- function(ls, p_value, ylimit) {
    sp_with_TF <- ls$sp_with_TF
    nonsp_with_TF <- ls$nonsp_with_TF
    sp_peaks_UDHS <- ls$sp_peaks_UDHS
    non_sp_peaks_UDHS <- ls$non_sp_peaks_UDHS

    a <- length(sp_with_TF)
    b <- length(sp_peaks_UDHS) - length(sp_with_TF)
    c <- length(nonsp_with_TF)
    d <- length(non_sp_peaks_UDHS) - length(nonsp_with_TF)

    plot_df <- data.frame(PeakType = rep(c("Spatial", "Non-Spatial"), each = 2),
                        KLF4 = rep(c("KLF4+", "KLF4-"), 2),
                        Count = c(a, b, c, d), 
                        Percentage = c(a/length(sp_peaks_UDHS)*100, b/length(sp_peaks_UDHS)*100, c/length(non_sp_peaks_UDHS)*100, d/length(non_sp_peaks_UDHS)*100))

    p_val <- p_value

    p <- ggplot(plot_df, aes(x = PeakType, y = Percentage, fill = PeakType)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), position = position_stack(vjust = 0.5), size = 4, color = "black", family = "Arial") +
    labs(title = paste0("p ", p_val), y = "Percentage", x = "Region Type") +
    scale_fill_manual(values = c("Spatial" = "#606060", "Non-Spatial" = "#E0E0E0")) +
    theme_bw() +
    ylim(0, ylimit) +
    theme(plot.title = element_text(size = 14, face = "bold", family = "Arial", hjust = 0.5), 
          legend.position = "none", 
          axis.title.x = element_text(size = 12, face = "bold", family = "Arial"), 
          axis.title.y = element_text(size = 12, face = "bold", family = "Arial"),
          axis.text = element_text(size = 10, family = "Arial"), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank())
}











##### Motif accessibility
motif_acc <- function(seu_obj, TF_name) {
    bed_path <- paste0("/project/zanglab_project/jw4xtu/motif_sites/mappable_results/", TF_name, "_map.bed")
    bed_files <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(bed_files) <- c("chr", "start", "end", "seq", "pval", "strand")
    motif_sites <- GRanges(seqnames = bed_files$chr, ranges = IRanges(start = bed_files$start + 1, end = bed_files$end))

    expression_matrix <- seu_obj@assays$peaks@counts
    peak <- rownames(expression_matrix)
    peak_info <- do.call(rbind, strsplit(peak, "-"))
    output_df <- data.frame(V1 = peak_info[,1], V2 = as.integer(peak_info[,2]), V3 = as.integer(peak_info[,3]))
    rownames(output_df) <- paste0(output_df$V1, "-", output_df$V2, "-", output_df$V3)
    input_peaks_gr <- GRanges(seqnames = output_df$V1, ranges = IRanges(start = output_df$V2, end = output_df$V3))

    overlaps <- findOverlaps(motif_sites, input_peaks_gr)
    input_peaks_overlapping <- input_peaks_gr[subjectHits(overlaps)]
    overlap_peaks <- paste0(as.character(seqnames(input_peaks_overlapping)), "-", start(input_peaks_overlapping), "-", end(input_peaks_overlapping))

    input_peaks <- output_df %>% mutate(region = rownames(output_df)) %>% dplyr::select(region)
    input_peaks <- input_peaks %>% mutate(rank = seq_len(nrow(.)))

    overlap_peak_rank <- input_peaks[input_peaks$region %in% overlap_peaks, ]

    peak_names <- overlap_peaks
    peak_matrix <- GetAssayData(seu_obj, assay = "peaks", slot = "counts")[peak_names, ]
    signal_per_spot <- Matrix::colMeans(peak_matrix)
    seu_obj$motif_accessibility <- signal_per_spot[colnames(seu_obj)]

    upper_lim <- quantile(signal_per_spot, 1, na.rm = TRUE)

    SpatialFeaturePlot(seu_obj, features = "motif_accessibility", pt.size.factor = 50) + 
        scale_fill_gradientn(colors = c("#FFFFFF", "#FDE5E5", "#FF3333", "#FF0000"), 
                             limits = c(0, upper_lim), oob = scales::squish) + 
        labs(fill = paste0(as.character(TF_name), " motif-containing\npeak accessibility")) + 
        theme(legend.position = "top", 
              legend.justification = c(0.5, 0),
              legend.box.spacing = unit(0, "pt"),
              legend.margin = margin(2,2,2,2),
              legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 12, family = "Arial"), 
              legend.key.width = unit(0.8, "cm"), 
              legend.key.height = unit(0.5, "cm"),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title = element_blank(),
              panel.grid = element_blank(), 
              panel.border = element_blank())
}
















##### Pearson's correlation between two TF list scaled ranks
correlation <- function(dt1, dt2, x, y) {
    dt1$scaled_rank <- round((max(dt1$rank) - dt1$rank) / (max(dt1$rank) - min(dt1$rank)), 3)
    dt2$scaled_rank <- round((max(dt2$rank) - dt2$rank) / (max(dt2$rank) - min(dt2$rank)), 3)

    sig_dt <- merge(dt1, dt2, by = "TF") %>% dplyr::select(TF, rank_avg_z_p_a_irwinhall_pvalue.x, scaled_rank.x, rank_avg_z_p_a_irwinhall_pvalue.y, scaled_rank.y)
    colnames(sig_dt) <- c("TF", "gene_pval", "gene_scaled_rank", "peak_pval", "peak_scaled_rank")

    dt <- sig_dt
    colnames(dt) <- c("TF", "pval1", "scale_rank_1", "pval2", "scale_rank_2")
    dt <- dt %>% mutate(overlap_TF = ifelse(pval1 < 0.05 & pval2 < 0.05, "True", "False"))

    correlation_res <- cor.test(dt$scale_rank_1, dt$scale_rank_2, method = "spearman")
    correlation <- round(correlation_res$estimate, 3)
    print(correlation)
    
    ggplot(dt, aes(x = scale_rank_1, y = scale_rank_2, color = overlap_TF)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = c("grey70", "#E41A1C")) + 
        theme_bw() + 
        theme(panel.grid.major = element_line(color = "grey85", linetype = "dotted"), 
              panel.grid.minor = element_blank(), 
              legend.position = "none",
              axis.text = element_text(size = 14, family = "Arial"), 
              axis.title = element_text(size = 16, face = "bold", family = "Arial")) + 
        labs(x = x, y = y, color = "TF Overlap")
}








##### Average accessiblity of peaksets
peak_avg_acc <- function(expression_matrix, peak, cell_df, title) {
    moran_expr <- Matrix::colMeans(expression_matrix[peak, , drop = FALSE])

    cell_df$moran_peak_avg_acc <- moran_expr

    ggplot(cell_df, aes(x = y, y = x, color = moran_peak_avg_acc)) +
        geom_point(size = 1.3) +
        scale_color_viridis_c(option = "D") +
        scale_y_reverse() +
        theme_bw() +
        labs(title = title, color = "Accessibility") +
        theme(legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
              legend.text = element_text(size = 12, family = "Arial"), 
              plot.title = element_text(size = 14, face = "bold", family = "Arial", hjust = 0.5),
              axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.y = element_blank(), 
              panel.grid.minor.y = element_blank())
}













##### Average accessibility of SVFs (similar as the function above)
SVF_avg_acc <- function(method, cell_df, data) {
  moran_expr <- Matrix::colMeans(obj$expression_matrix[data, , drop = FALSE])
  cell_df$avg_acc <- moran_expr
  cell_df$method <- method
  return(cell_df)
}

SVF_avg_acc_plot <- function(cell_df, alternative) {
    cell_df$cell_type <- factor(cell_df$cell_type, levels = c("Radial glia", "Postmitotic premature neurons"))

    plot <- ggplot(cell_df, aes(x = method, y = avg_acc, fill = cell_type)) +
            geom_boxplot(position = position_dodge(width = 0.7), width = 0.6) +
            scale_fill_manual(values = c("Radial glia" = "#E64B35", "Postmitotic premature neurons" = "gray70")) +
            labs(x = "", y = "Average accessibility", fill = "Cell type") +
            theme_bw() +
            theme(legend.position = "right", 
                  legend.title = element_text(size = 14, face = "bold", family = "Arial"), 
                  legend.text = element_text(size = 12, family = "Arial"), 
                  axis.text = element_text(size = 14, family = "Arial"), 
                  axis.title = element_text(size = 16, face = "bold", family = "Arial"))

    # Wilcoxon per method
    annotation_label <- function(p) {
    if (p < 0.00001) "*****"
    else if (p < 0.0001) "****"
    else if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else "NS"
    }

    stats <- cell_df %>%
                group_by(method) %>%
                summarize(p = wilcox.test(avg_acc ~ cell_type, exact = FALSE, alternative = alternative)$p.value,
                          y = max(avg_acc, na.rm = TRUE) - 0.01, .groups = "drop") %>%
                mutate(star = vapply(p, annotation_label, character(1)))

    p_stat <- plot + geom_text(data = stats, aes(x = method, y = y, label = star), vjust = 0, size = 3, inherit.aes = FALSE)

    return(p_stat)
}