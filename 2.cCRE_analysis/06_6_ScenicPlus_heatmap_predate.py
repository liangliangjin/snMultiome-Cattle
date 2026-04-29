# Process scplusmdata.h5mu (output from new ScenicPlus snakemake pipeline) for heatmap_dotplot
import re
import numpy as np
import pandas as pd
import mudata as md
from scenicplus.RSS import regulon_specificity_scores

h5mu_path = "outs/scplusmdata.h5mu"
group_key = "main"
auc_modality = "direct_gene_based_AUC"
rna_modality = "scRNA_counts"
metadata_key = "direct_e_regulon_metadata"
out_csv = "scenicplus_TF_expression_AUC_RSS_long_table.csv"

mdata = md.read_h5mu(h5mu_path)

rna = mdata[rna_modality]
auc = mdata[auc_modality]
common_cells = rna.obs_names.intersection(auc.obs_names)
rna = rna[common_cells].copy()
auc = auc[common_cells].copy()
groups = rna.obs[group_key].astype(str)

auc_df = auc.to_df()
mean_auc = auc_df.groupby(groups).mean().T
expr_df = rna.to_df()  # cells x genes
mean_expr = expr_df.groupby(groups).mean().T

meta = mdata.uns[metadata_key].copy()
regulon_info = (
    meta[["Gene_signature_name", "TF", "regulation"]]
    .dropna()
    .drop_duplicates()
    .groupby("Gene_signature_name")
    .agg({
        "TF": "first",
        "regulation": "first"
    })
)
regulon_to_tf = regulon_info["TF"]
regulon_to_regulation = regulon_info["regulation"]

common_regulons = mean_auc.index.intersection(regulon_to_tf.index)
common_regulons = [
    r for r in common_regulons
    if regulon_to_tf.loc[r] in mean_expr.index
]

mean_auc = mean_auc.loc[common_regulons]
regulon_to_tf = regulon_to_tf.loc[common_regulons]
regulon_to_regulation = regulon_to_regulation.loc[common_regulons]
tf_expr_for_regulon = pd.DataFrame(
    index=mean_auc.index,
    columns=mean_auc.columns,
    dtype=float
)

for regulon in mean_auc.index:
    tf = regulon_to_tf.loc[regulon]
    tf_expr_for_regulon.loc[regulon] = mean_expr.loc[tf, mean_auc.columns]


# Calculate RSS
rss = regulon_specificity_scores(
    scplus_mudata=mdata,
    variable=f"{rna_modality}:{group_key}",
    modalities=[auc_modality]
)

if isinstance(rss, dict):
    rss = rss[auc_modality]

rss = rss.loc[mean_auc.columns, mean_auc.index].T

# Scale AUC, TF expression, RSS
def row_zscore(df):
    df = df.astype(float)
    out = df.sub(df.mean(axis=1), axis=0)
    out = out.div(df.std(axis=1).replace(0, np.nan), axis=0)
    return out.fillna(0)

def row_minmax(df):
    df = df.astype(float)
    out = df.sub(df.min(axis=1), axis=0)
    out = out.div((df.max(axis=1) - df.min(axis=1)).replace(0, np.nan), axis=0)
    return out.fillna(0)

auc_scale = row_zscore(mean_auc).clip(-2.5, 2.5)
expr_scale = row_zscore(tf_expr_for_regulon).clip(-2.5, 2.5)
rss_scale = row_minmax(rss)

# Build long table
records = []
for regulon in mean_auc.index:
    tf = regulon_to_tf.loc[regulon]
    regulation = regulon_to_regulation.loc[regulon]
    for cell_type in mean_auc.columns:
        records.append({
            "Regulon": regulon,
            "TF": tf,
            "Cell_type": cell_type,
            "mean_AUC": mean_auc.loc[regulon, cell_type],
            "AUC_scale": auc_scale.loc[regulon, cell_type],
            "Expression": tf_expr_for_regulon.loc[regulon, cell_type],
            "Expression_scale": expr_scale.loc[regulon, cell_type],
            "RSS": rss.loc[regulon, cell_type],
            "RSS_scale": rss_scale.loc[regulon, cell_type],
            "regulation": regulation,
            "Repressor_activator": "Activators" if regulation == 1 else "Repressors"
        })

rel_data = pd.DataFrame(records)
print(regulon_to_regulation.value_counts())
print(pd.crosstab(rel_data["regulation"], rel_data["Repressor_activator"]))
# order regulons by max RSS cell type
rss_order = []
cell_type_order = list(mean_auc.columns)

for cell_type in cell_type_order:
    regs = (
        rel_data[rel_data["Cell_type"] == cell_type]
        .sort_values("RSS_scale", ascending=False)["Regulon"]
        .tolist()
    )
    for r in regs:
        if r not in rss_order:
            rss_order.append(r)

rel_data["Regulon"] = pd.Categorical(
    rel_data["Regulon"],
    categories=rss_order,
    ordered=True
)

rel_data["Cell_type"] = pd.Categorical(
    rel_data["Cell_type"],
    categories=cell_type_order,
    ordered=True
)

rel_data.to_csv(out_csv, index=False)