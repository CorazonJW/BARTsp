'''
conda create -n scripro python=3.8
conda activate scripro
conda install -c liulab-dfci lisa2
lisa install hg38 oneshot ./hg38_scripro.h5 --force
lisa install mm10 oneshot ./mm10_scripro.h5 --force
pip install scripro
scripro install_reference -i TF_target_RP.h5
'''


##### 1. Prepare input
import scanpy as sc

# Convert to h5ad file
# Expression matrix
adata = sc.read_text("~/expression_matrix.txt")
adata_all = adata.transpose().copy()
adata_all.obs_names = adata_all.obs_names.str.replace(r"\.", "-", regex=True)
# Metadata
df_annot = pd.read_csv("~/cell_metadata.csv", index_col=0)
common_barcodes = adata_all.obs_names.intersection(df_annot.index)
adata_all.obs.loc[common_barcodes, "cell_type"] = df_annot.loc[common_barcodes, "cell_type"]
adata_all.obs["cluster"] = df_annot["cell_type"].values
adata_sub = adata_all[adata_all.obs["cluster"].isin(["Enterocyte_Progenitor", "Enterocyte_Immature", "Enterocyte_Mature"])].copy()

adata_sub.write("~/mm_small_intestine_expression_matrix.h5ad", compression="gzip")

# Check
adata_check = sc.read_h5ad("~/mm_small_intestine_expression_matrix.h5ad")
print(adata_check)



##### 2. Run SCRIPro
'''
scripro enrich_rna -i mm_small_intestine_expression_matrix.h5ad -n 50 -s mm10 -p mm_Visium_HD -t 8
'''



##### 3. Post-processing SCRIPro results
# Order cells based on pseudotime
import pandas as pd
import numpy as np

df = pd.read_pickle("~/scripro_result.pkl")
tf_activity_matrix = df.tf_score
supercell_info = df.Ori_Data.cellgroup

# Annotate SuperCell with cell type label
df_annot = pd.read_csv("~/cell_metadata.csv", index_col=0)
df_supercell = supercell_info.join(df_annot[['cell_type']], how='left')

# Purity check: SuperCell identity defined by the major cell type
ct_tab = (
    df_supercell.groupby(['new_leiden', 'cell_type'])
      .size()
      .unstack(fill_value=0)
)

majority_label = ct_tab.idxmax(axis=1)
purity = ct_tab.max(axis=1) / ct_tab.sum(axis=1)

summary = pd.DataFrame({
    'n_cells': ct_tab.sum(axis=1),
    'majority_cell_type': majority_label,
    'purity': purity
}).sort_values('n_cells', ascending=False)

summary['is_mixed'] = summary['purity'] < 0.75
summary_sorted = summary.sort_values(by='purity', ascending=False)
print(summary_sorted)

# Summarize TF activity by cell type
sc2ct = summary_sorted['majority_cell_type'].to_dict()
tmp = tf_activity_matrix.copy()
tmp['cell_type'] = tmp.index.to_series().map(sc2ct).fillna('unknown')
tf_activity_matrix = (tmp.groupby('cell_type').mean())

# Obtain p-value matrix and keep the smallest for one TF at each stage
rows = []
for sc, res in df.results.items():
    tmp = res[['factor', 'summary_p_value']].copy()
    tmp.columns = ['TF', 'pval']
    tmp['supercell'] = sc
    rows.append(tmp)

pval_long = pd.concat(rows, ignore_index=True)
pval_matrix = pval_long.pivot_table(index='supercell', columns='TF', values='pval', aggfunc='min')

tmp_p = pval_matrix.copy()
tmp_p['cell_type'] = tmp_p.index.to_series().map(sc2ct).fillna('unknown')
pval_matrix = (tmp_p.groupby('cell_type').min())

# Divide TF into upstream and downstream-active
upstream_act = tf_activity_matrix.loc['Enterocyte_Immature']
downstream_act = tf_activity_matrix.loc['Enterocyte_Mature']
upstream_p = pval_matrix.loc['Enterocyte_Immature']
downstream_p = pval_matrix.loc['Enterocyte_Mature']

common_tfs_up = upstream_act.index.intersection(upstream_p.index)
upstream_result = pd.DataFrame({'upstream_activity': upstream_act.loc[common_tfs_up], 'p-value': upstream_p.loc[common_tfs_up], 'label': 'upstream'})
common_tfs_down = downstream_act.index.intersection(downstream_p.index)
downstream_result = pd.DataFrame({'downstream_activity': downstream_act.loc[common_tfs_down], 'p-value': downstream_p.loc[common_tfs_down], 'label': 'downstream'})

upstream_result.to_csv("~/scripro_results_upstream.csv")
downstream_result.to_csv("~/scripro_results_downstream.csv")


