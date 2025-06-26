# conda activate SpaTrack

import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import spaTrack as spt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sc.settings.verbosity = 0
plt.rcParams['figure.dpi'] = 200

##### 1. Prepare inputs
# Expression matrix
adata = sc.read_text("~/data/expression_matrix.txt")
adata_all = adata.transpose().copy()
adata_all.obs_names = adata_all.obs_names.str.replace(r"\.", "-", regex=True)

# Metadata  
df_annot = pd.read_csv("~/data/cell_metadata.csv", index_col=0)
common_barcodes = adata_all.obs_names.intersection(df_annot.index)
adata_all.obs.loc[common_barcodes, "cell_type"] = df_annot.loc[common_barcodes, "cell_type"]

# Spatial coordinates
spatial_coords = pd.read_csv("~/data/spatial_coordinates.csv", index_col=0)
spatial_coords = spatial_coords.reset_index().rename(columns={"index": "barcode"})
spatial_coords = spatial_coords.set_index("barcode")
spatial_coords = spatial_coords.loc[adata_all.obs_names]

# Create adata
adata_all.obs["cluster"] = df_annot["cell_type"].values
adata_all.obsm["X_spatial"] = spatial_coords[["x", "y"]].values
adata_all.obs[["x", "y"]] = adata_all.obsm["X_spatial"]
adata_all.layers["counts"] = adata_all.X
adata_all.obs['CellID']=adata_all.obs.index
adata_sub = adata_all[adata_all.obs["cluster"].isin(["core", "edge", "transitory"])].copy()

sc.pp.filter_genes(adata_sub, min_cells=10)
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)
sc.pp.calculate_qc_metrics(adata_sub, percent_top=None, log1p=False, inplace=True)

##### 2. Construct trajectory
start_cells = spt.set_start_cells(adata_sub, select_way='cell_type', cell_type='core')
adata_sub.obsp["trans"] = spt.get_ot_matrix(adata_sub, data_type="spatial", alpha1=0.6, alpha2=0.4)
adata_sub.obs["ptime"] = spt.get_ptime(adata_sub, start_cells)
adata_sub.uns["E_grid"], adata_sub.uns["V_grid"] = spt.get_velocity(adata_sub, basis="spatial", n_neigh_pos=50)

# Plot
fig, axs = plt.subplots(figsize=(5, 5))
sc.pl.embedding(adata_sub, basis='spatial', color='ptime', show=False, ax=axs, color_map='Reds', title='ptime', size=100)
if len(fig.axes) > 1:
    fig.delaxes(fig.axes[1])
axs.xaxis.set_major_locator(ticker.MultipleLocator(500))
axs.yaxis.set_major_locator(ticker.MultipleLocator(500))
plt.savefig("~/results/spaTrack_ptime_spatial_plot.png", dpi=300, bbox_inches='tight')


##### 3. Predict trajectory-related genes
sub_adata_path = adata_sub[adata_sub.obs['cluster'].isin(["core", "edge", "transitory"])]
sub_adata_path = sub_adata_path[:, ~sub_adata_path.var_names.duplicated()].copy()
sub_adata_path = sub_adata_path[sub_adata_path.obs["ptime"].argsort()].copy()
sub_adata_path = spt.filter_gene(sub_adata_path, min_exp_prop=0.1, hvg_gene=3000)

df_res  = spt.ptime_gene_GAM(sub_adata_path, core_number=8)
df_sig_res = df_res.loc[(df_res['model_fit']>0.05) & (df_res['fdr']<0.05) & (df_res['pvalue']<0.05)]
df_sig_res.to_csv("~/results/spatrack_sort_exp_sig_with_trend.tsv", sep="\t")

sort_exp_sig = spt.order_trajectory_genes(sub_adata_path, df_sig_res, cell_number=20)
spt.plot_trajectory_gene_heatmap(sort_exp_sig,smooth_length=100,gene_label_size=15,cmap_name='twilight_shifted')
sort_exp_sig.to_csv("~/results/spatrack_sort_exp_sig.tsv", sep="\t")



##### 4. Predict trajectory-related transcription factors
adata_sub = adata_sub[:, ~adata_sub.var_names.duplicated()].copy()
adata_sub.obs["cell_type"] = adata_sub.obs["cell_type"].astype(str)
adata_sub.obs = adata_sub.obs.astype(str)
adata_sub.write("~/results/spatrack_dedup.adata.h5ad")

ptime_df = adata_sub.obs[["ptime"]]
ptime_df.to_csv("~/results/pseudotime.tsv", sep="\t")

# only use genes with decreased expression to predict upstream-active TFs
adata = sc.read_h5ad("~/results/spatrack_dedup.adata.h5ad")
DEGs = pd.read_csv("~/results/spatrack_sort_exp_sig_with_trend.tsv", sep = "\t")
DEGs_decreasing = DEGs[DEGs['pattern'] == 'decrease']
decreasing_genes = DEGs_decreasing["gene"].values
common_genes = adata.var_names.intersection(decreasing_genes)
adata_decreasing = adata[:, common_genes].copy()
adata_decreasing.write("~/results/spatrack_decrease.adata.h5ad")

# run TF prediction model
gr = spt.Trainer(data_type="p_time", expression_matrix_path="~/data/spatrack_decrease.adata.h5ad", tfs_path="~/data/hs.TF.txt", ptime_path="~/results/pseudotime.tsv", min_cells=0)
import numpy as np; gr.output_data = np.nan_to_num(gr.output_data)
print(gr.input_data.shape, gr.output_data.shape)
print(gr.output_data.std(axis=0).mean())
gr.get_dataloader()

gr.run()
gr.network_df
gr.network_df.to_csv('~/results/spatrack_TF_prediction_spatrack_TF_decrease.csv')


##### 5. Check SpaTrack prediction performance
'''
# Continue in R
# Check whether identify known functional TFs
spatrack_result <- read.csv("~/results/spatrack_TF_prediction_spatrack_TF_decrease.csv")

spatrack_result$rank <- seq_len(nrow(spatrack_result))

sig_result <- spatrack_result %>% filter(weight > cutoff)
tf_interest <- c("GRHL3", "TCF4", 
                 "TP63", "GRHL2", "SOX2", "KLF4", "TP73", 
                 "SNAI2", "ZEB1",  "FOSL1", "STAT3", 
                 "FOS", "JUN", "JUND", "TP53")

sig_TF <- sig_result %>% filter(TF %in% tf_interest)


# Check the rank difference between functional and other TFs
library(stringr)
library(ggsignif)
library(ggplot2)
spatrack_result$TF_type <- ifelse(spatrack_result$TF %in% tf_interest, "functional", "not_functional")
spatrack_result$TF_type <- factor(spatrack_result$TF_type, levels = c("functional", "not_functional"))
spatrack_result$rank <- as.numeric(spatrack_result$rank)
spatrack_result$scaled_rank <- round((max(spatrack_result$rank) - spatrack_result$rank) / (max(spatrack_result$rank) - min(spatrack_result$rank)), 3)

wilcox_result <- wilcox.test(spatrack_result$scaled_rank ~ spatrack_result$TF_type, data = spatrack_result, exact = FALSE)
p_val <- round(wilcox_result$p.value, 2)

if (p_val < 0.1) {
  annotation_label <- if (p_val < 0.01) "***"
                    else if (p_val < 0.05) "**"
                    else "*"
} else {
  annotation_label <- "NS"
}

p <- ggplot(spatrack_result, aes(x = TF_type, y = scaled_rank, fill = TF_type)) +
      geom_boxplot() +
      geom_signif(comparisons = list(c("functional", "not_functional")), annotations = annotation_label, map_signif_level = TRUE, textsize = 5, vjust = 1.5) +
      theme_bw() +
      labs(title = "", x = "TF type", y = "Scaled rank") +
      theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 12))
ggsave("~/results/spatrack_TF_rank_spatrack_TF.pdf", p, width = 3.5, height = 3.5)
'''