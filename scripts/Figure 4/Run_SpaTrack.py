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
adata = sc.read_text("~/exprMatrix.tsv.gz")
adata_all = adata.transpose().copy()
adata_all.var_names = adata_all.var_names.str.replace('^atac-', '', regex=True)
# Metadata
df_annot = pd.read_csv("~/E13_processed_metadata.tsv", sep="\t")
common_barcodes = adata_all.obs_names.intersection(df_annot.index)
adata_all.obs.loc[common_barcodes, "predicted.id"] = df_annot.loc[common_barcodes, "predicted.id"]
# Spatial coordinates
spatial_coords = pd.read_csv("~/Spots.coords.tsv", sep="\t", header=None, names=["barcode", "x", "y"])
spatial_coords = spatial_coords.set_index("barcode")
spatial_coords = spatial_coords.loc[adata_all.obs_names]

# Create adata
adata_all.obs["cluster"] = df_annot["predicted.id"].values
adata_all.obsm["X_spatial"] = spatial_coords[["x", "y"]].values
adata_all.obs[["x", "y"]] = adata_all.obsm["X_spatial"]
adata_all.layers["counts"] = adata_all.X
adata_all.obs['CellID']=adata_all.obs.index
adata_sub = adata_all[adata_all.obs["cluster"].isin(["Radial glia", "Postmitotic premature neurons"])].copy()

sc.pp.filter_genes(adata_sub, min_cells=10)
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)
sc.pp.calculate_qc_metrics(adata_sub, percent_top=None, log1p=False, inplace=True)



##### 2. Construct trajectory
start_cells = spt.set_start_cells(adata_sub, select_way='cell_type', cell_type='Radial glia')
adata_sub.obsp["trans"] = spt.get_ot_matrix(adata_sub, data_type="spatial", alpha1=0.6, alpha2=0.4)
adata_sub.obs["ptime"] = spt.get_ptime(adata_sub, start_cells)
adata_sub.uns["E_grid"], adata_sub.uns["V_grid"] = spt.get_velocity(adata_sub, basis="spatial", n_neigh_pos=50)

# Plot
fig, axs = plt.subplots(figsize=(5, 5))
sc.pl.embedding(adata_sub, basis='spatial', color='ptime', show=False, ax=axs, color_map='Reds', title='ptime', size=100)
axs.xaxis.set_major_locator(ticker.MultipleLocator(500))
axs.yaxis.set_major_locator(ticker.MultipleLocator(500))
plt.savefig("~/ptime_spatial_plot.png", dpi=300, bbox_inches='tight')



##### 3. Predict trajectory-related genes
sub_adata_path = adata_sub[adata_sub.obs['cluster'].isin(["Radial glia", "Postmitotic premature neurons"])]
sub_adata_path = sub_adata_path[:, ~sub_adata_path.var_names.duplicated()].copy()
sub_adata_path = sub_adata_path[sub_adata_path.obs["ptime"].argsort()].copy()
sub_adata_path = spt.filter_gene(sub_adata_path, min_exp_prop=0.1, hvg_gene=3000)

df_res  = spt.ptime_gene_GAM(sub_adata_path, core_number=8)
df_sig_res = df_res.loc[(df_res['model_fit']>0.6) & (df_res['fdr']<0.05) & (df_res['pvalue']<0.05)]
df_sig_res.to_csv("~/spatrack_sort_exp_sig_with_trend.tsv", sep="\t")

sort_exp_sig = spt.order_trajectory_genes(sub_adata_path, df_sig_res, cell_number=20)
spt.plot_trajectory_gene_heatmap(sort_exp_sig,smooth_length=100,gene_label_size=15,cmap_name='twilight_shifted')
sort_exp_sig.to_csv("~/spatrack_sort_exp_sig.tsv", sep="\t")



##### 4. Predict trajectory-related transcription factors
ptime_df = adata_sub.obs[["ptime"]]
ptime_df.to_csv("~/pseudotime.tsv", sep="\t")

adata_sub = adata_sub[:, ~adata_sub.var_names.duplicated()].copy()
adata_sub.write("~/spatrack_dedup.adata.h5ad")

# only use genes with decreased expression to predict upstream-active TFs
adata = sc.read_h5ad("~/spatrack_dedup.adata.h5ad")
DEGs = pd.read_csv("~/spatrack_sort_exp_sig_with_trend.tsv", sep = "\t")
DEGs_decreasing = DEGs[DEGs['pattern'] == 'decrease']
decreasing_genes = DEGs_decreasing["gene"].values
common_genes = adata.var_names.intersection(decreasing_genes)
adata_decreasing = adata[:, common_genes].copy()
adata_decreasing.write("~/spatrack_decrease.adata.h5ad")

# only use genes with increased expression to predict upstream-active TFs
DEGs_increasing = DEGs[DEGs['pattern'] == 'increase']
increasing_genes = DEGs_increasing["gene"].values
common_genes = adata.var_names.intersection(increasing_genes)
adata_increasing = adata[:, common_genes].copy()
adata_increasing.write("~/spatrack_increase.adata.h5ad")

gr = spt.Trainer(data_type="p_time", expression_matrix_path="~/spatrack_decrease.adata.h5ad", tfs_path="~/SpaTrack_ms.TF.tsv", ptime_path="~/pseudotime.tsv", min_cells=0)
gr = spt.Trainer(data_type="p_time", expression_matrix_path="~/spatrack_increase.adata.h5ad", tfs_path="~/SpaTrack_ms.TF.tsv", ptime_path="~/pseudotime.tsv", min_cells=0)

import numpy as np; gr.output_data = np.nan_to_num(gr.output_data)
print(gr.input_data.shape, gr.output_data.shape)
print(gr.output_data.std(axis=0).mean())
gr.get_dataloader()

gr.run()
gr.network_df
gr.network_df.to_csv('~/spatrack_results_upstream.csv')
gr.network_df.to_csv('~/spatrack_results_downstream.csv')

