# Using pycisTopic on cattle multiome data
# Some of the steps differ from the tutorial because certain information has already been generated through ArchR/Seurat.
import os
from scipy.io import mmread
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import numpy as np
import pickle
import resource

# 0. Setting up the analysis directories & common variables
# paths
work_path = "/storage/public/home/2021060195/03.snATAC/06.scenicplus/"
#os.makedirs(out_path, exist_ok = True)
fragments_path = "/storage/public/home/2021060195/SC/combine_all/"
out_dir = "outs"
#os.makedirs(out_dir, exist_ok = True)

# 1. Reading data & creating adata object
# Read metadata
metadata_file = os.path.join(work_path, 'cell_metadata.tsv')
metadata = pd.read_csv(metadata_file, sep = '\t', index_col = 0)
print(metadata.head())
print(metadata.shape)

# Load the ATAC sparse matrix
atac_path = os.path.join(work_path, "atac_counts.mtx")
atac_matrix = mmread(atac_path).tocsr()
genes_path = os.path.join(work_path, "peaks_regions.tsv")
cells_path = os.path.join(work_path, "cell_names.tsv")
genes = pd.read_csv(genes_path, header=None, sep='\t')[0].values
cells = pd.read_csv(cells_path, header=None, sep='\t')[0].values
atac_counts_df = pd.DataFrame.sparse.from_spmatrix(atac_matrix, index=genes, columns=cells)
print(atac_counts_df.head())

# Load the RNA sparse matrix
rna_path = os.path.join(work_path, "rna_counts.mtx")
rna_matrix = mmread(rna_path).tocsr()
genes_path = os.path.join(work_path, "gene_name.tsv")
cells_path = os.path.join(work_path, "cell_names.tsv")
genes = pd.read_csv(genes_path, header=None, sep='\t')[0].values
cells = pd.read_csv(cells_path, header=None, sep='\t')[0].values
rna_counts_df = pd.DataFrame.sparse.from_spmatrix(rna_matrix, index=genes, columns=cells)
print(rna_counts_df.head())

# creating adata object
adata = sc.AnnData(X=rna_counts_df.T, obs=metadata, var=pd.DataFrame(index=rna_counts_df.index))
print(adata)

columns_to_keep = ['sample', 'cattle', 'main', 'tissue', 'celltype_cattle', 'WNN_Clusters']
adata.obs = adata.obs[columns_to_keep]

adata.obs['barcode'] = adata.obs.index.str.split('#').str[1]
adata.obs.index = adata.obs['barcode'] + '-' + adata.obs['sample'] + '___' + adata.obs['sample']

adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(adata)
# save
adata.write("adata.h5ad")

# 2. pseudobulk generation
fragments_dict = {sample: os.path.join(fragments_path, f"{sample}.gz") for sample in adata.obs['sample'].unique()}
print(fragments_dict)

columns = ['chrom', 'start', 'end', 'barcode', 'count']
all_fragments_df = pd.DataFrame(columns=columns)
sample_to_cells = metadata.groupby('sample').apply(lambda x: x.index.tolist()).to_dict()
print(sample_to_cells)

fragments_df = pd.read_csv("fragments_paths.csv", index_col=0)
fragments_dict = fragments_df['inputFiles'].to_dict()

cell_data = pd.read_table("cell_metadata.tsv", index_col = 0)
cell_data['barcode'] = cell_data.index.str.split('#').str[1]
cell_data.index = cell_data['barcode'] + '-' + cell_data['Sample']
cell_data.head()

chromsizes = pd.read_table("ARSUCD1.2.chr.sizes", header = None, names = ["Chromosome", "End"])
chromsizes.insert(1, "Start", 0)
chromsizes.head()

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

bw_paths, bed_paths = export_pseudobulk(
   input_data = cell_data,
   variable = "main",
   sample_id_col = "Sample",
   chromsizes = chromsizes,
   bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
   bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
   path_to_fragments = fragments_dict,
   n_cpu = 10,
   normalize_bigwig = True,
   temp_dir = "/tmp",
   split_pattern = "-"
)

with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
   for v in bw_paths:
       _ = f.write(f"{v}\t{bw_paths[v]}\n")
        
        
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
   for v in bed_paths:
       _ = f.write(f"{v}\t{bed_paths[v]}\n")
        

# 3. Creating a cisTopic object
path_to_regions = os.path.join("peaks_edgeR.bed") #Or select the regions of interest.
path_to_blacklist = os.path.join("ARSUCD1.2_blacklist.bed")
pycistopic_qc_output_dir = "outs/qc"

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)
    

cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)

cistopic_obj.merge(cistopic_obj_list[1:])
print(cistopic_obj)

import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_celltype_obj.pkl"), "wb")
)

# Adding cell information
import os
import pickle
with open("adata.pkl", "rb") as file:
    adata = pickle.load(file)
	

with open(os.path.join(out_dir, "cistopic_celltype_obj.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)
##Keep the barcode format consistent.
print(cistopic_obj.cell_data.head())
print(adata.obs.head())

cistopic_obj.add_cell_data(adata.obs)
cells_with_adata = cistopic_obj.cell_data.index[
    cistopic_obj.cell_data['main'].notna()
]
len(cells_with_adata) == len(set(adata.obs.index))
cistopic_obj.subset(cells = cells_with_adata)
len(cistopic_obj.cell_data)
# save new cistopic_obj
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_celltype_obj.pkl"), "wb")
)


# 4. Run models
with open(os.path.join(out_dir, "cistopic_celltype_obj.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)

os.environ['MALLET_MEMORY'] = '900G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="Mallet-202108/bin/mallet"
#os.mkdir(out_dir+'mallet', exist_ok=True)
# Test the appropriate number of topics.
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,105,110,115,120,125,130,135,140,145,150],
    n_cpu=60,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(out_dir, "mallet"),
    save_path=None,
    mallet_path=mallet_path,
)
	
model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + '/model_selection.pdf')
					 
# Set the n_topics to 75.
os.environ['MALLET_MEMORY'] = '900G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[80],
    n_cpu=30,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(out_dir, "mallet"),
    save_path=os.path.join(out_dir, "mallet"),
    mallet_path=mallet_path,
)
	
	
pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

# # Add model to cisTopicObject
with open(os.path.join(out_dir, "models.pkl"), "rb") as file:
    models = pickle.load(file)

with open(os.path.join(out_dir, "cistopic_celltype_obj.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)
	

from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 75,
    return_model = True
)

cistopic_obj.add_LDA_model(model)
print(cistopic_obj.selected_model)

cistopic_obj.cell_data['main'] = cistopic_obj.cell_data['main'].astype(str)
# save
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_celltype_with_model.pkl"), "wb")
)

# 5. Clustering
# It is only used to compare the classification of the models. The subsequent cell types are based on the annotations in the previous ArchR/Seurat.
import pickle
with open(os.path.join(out_dir, "cistopic_celltype_with_model.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)

# cisTopic topics broadly recapitulated major cell groups
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [1],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)
# Check if projections attribute exists and is a dictionary
if not hasattr(cistopic_obj, 'projections') or not isinstance(cistopic_obj.projections, dict):
    cistopic_obj.projections = {}

# Check if 'cell' key exists in projections dictionary and initialize it if necessary
if 'cell' not in cistopic_obj.projections:
    cistopic_obj.projections['cell'] = {}
run_umap(
    cistopic_obj,
    target='cell',
    scale=True,
    n_neighbors=30,
    min_dist=0.5
)

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['main'],
    target='cell', num_columns=4,
    text_size=10,
    dot_size=5)
plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=5
)
#save
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_celltype_umap.pkl"), "wb")
)

# 6. Topic binarization & QC
from pycisTopic.topic_binarization import binarize_topics
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=False, num_columns=5
)
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=False, num_columns=5
)
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=False,
    num_columns=5, nbins=100)

from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img
topic_qc_metrics = compute_topic_metrics(cistopic_obj)
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
fig.savefig(out_dir + 'topic_binarization/Topic_qc.pdf', bbox_inches='tight')
plt.show()

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='main',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)
topic_annot
# Save
pickle.dump(
    topic_qc_metrics,
    open(os.path.join(out_dir, "topic_binarization/Topic_qc_metrics_annot.pkl"), "wb")
)
pickle.dump(
    region_bin_topics_top_3k,
    open(os.path.join(out_dir, "topic_binarization/region_bin_topics_top_3k.pkl"), "wb")
)
pickle.dump(
    region_bin_topics_otsu,
    open(os.path.join(out_dir, "topic_binarization/region_bin_topics_otsu.pkl"), "wb")
)
pickle.dump(
    binarized_cell_topic,
    open(os.path.join(out_dir, "topic_binarization/binarized_cell_topic.pkl"), "wb")

# 7. Differentially Accessible Regions (DARs)
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)
# save
pickle.dump(
    imputed_acc_obj,
    open(os.path.join(out_dir, "imputed_acc_obj.pkl"), "wb")
)

variable_regions = find_highly_variable_features(
    imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)
len(variable_regions)

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='main',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=16,
    _temp_dir="/storage/public/home/2021060195/temp",
    split_pattern = '-'
)
# save
pickle.dump(
    markers_dict,
    open(os.path.join(out_dir, "markers_dict.pkl"), "wb")
)

#
with open(os.path.join(out_dir, "markers_dict.pkl"), "rb") as file:
    markers_dict = pickle.load(file)
	
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")



# save three sets of genomic regions
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)
from pycisTopic.utils import region_names_to_coordinates
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )
	
# 8. Gene activity
# Load cisTopic object
import pickle
with open(os.path.join(out_dir, "cistopic_celltype_umap.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)
# Load imputed accessibility
with open(os.path.join(out_dir, "imputed_acc_obj.pkl"), "rb") as file:
    imputed_acc_obj = pickle.load(file)
# Load DARs
with open(os.path.join(out_dir, "markers_dict.pkl"), "rb") as file:
    markers_dict = pickle.load(file)
# Get TSS annotations
import pybiomart as pbm
import pyranges as pr
# For cattle
#dataset = pbm.Dataset(name='btaurus_gene_ensembl', host='http://www.ensembl.org')
annot = dataset.query(attributes=['chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name', 'transcription_start_site', 'transcript_biotype'])
annot.columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene','Transcription_Start_Site', 'Transcript_type']
annot['Chromosome'] = 'chr' + annot['Chromosome'].astype(str)
annot = annot[annot.Transcript_type == 'protein_coding']
annot['Strand'] = annot['Strand'].astype(str)
annot.loc[annot['Strand'] == '1', 'Strand'] = '+'
annot.loc[annot['Strand'] == '-1', 'Strand'] = '-'
pr_annotation = pr.PyRanges(annot.dropna(axis = 0))
pickle.dump(
    pr_annotation,
    open(os.path.join("pr_annotation.pkl"), "wb")
)

import requests
target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.chrom.sizes'
chromsizes = pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
chromsizes=pr.PyRanges(chromsizes)
pickle.dump(
    chromsizes,
    open(os.path.join("chromsizes.pkl"), "wb")
)


from pycisTopic.gene_activity import get_gene_activity
gene_act, weigths = get_gene_activity(
    imputed_acc_obj,
    pr_annotation,
    chromsizes,
    use_gene_boundaries=True,
    upstream=[1000, 100000],
    downstream=[1000,100000], 
    distance_weight=True,
    decay_rate=1,
    extend_gene_body_upstream=10000,
    extend_gene_body_downstream=500,
    gene_size_weight=False,
    gene_size_scale_factor='median',
    remove_promoters=False,
    average_scores=True,
    scale_factor=1,
    extend_tss=[10,10],
    gini_weight = True, 
    return_weights= True,
    project='Gene_activity')

pickle.dump(
    gene_act,
    open(os.path.join(out_dir, "Gene_activity.pkl"), "wb")
)

DAG_markers_dict= find_diff_features(
    cistopic_obj,
    gene_act,
    variable='main',
    var_features=None,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=20,
    _temp_dir="/storage/public/home/2021060195/temp",
    split_pattern = '-')

pickle.dump(
    DAG_markers_dict,
    open(os.path.join(out_dir, "DAG_markers_dict.pkl"), "wb")
)
print("Number of DAGs found:")
print("---------------------")
for x in DAG_markers_dict:
    print(f"  {x}: {len(DAG_markers_dict[x])}")

# 9. Exporting to loom
# Load cisTopic object
import pickle
with open(os.path.join(out_dir, "cistopic_celltype_umap.pkl"), "rb") as file:
    cistopic_obj = pickle.load(file)

# Load imputed accessibility
with open(os.path.join(out_dir, "imputed_acc_obj.pkl"), "rb") as file:
    imputed_acc_obj = pickle.load(file)

# Load region binarized topics
with open(os.path.join(out_dir, "topic_binarization/region_bin_topics_otsu.pkl"), "rb") as file:
    region_bin_topics_otsu = pickle.load(file)

# Load cell binarized topics
with open(os.path.join(out_dir, "topic_binarization/binarized_cell_topic.pkl"), "rb") as file:
    binarized_cell_topic = pickle.load(file)

# Load DARs
with open(os.path.join(out_dir, "markers_dict.pkl"), "rb") as file:
    markers_dict = pickle.load(file)

# Load DARs
with open(os.path.join(out_dir, "DAG_markers_dict.pkl"), "rb") as file:
    DAG_markers_dict = pickle.load(file)


from pycisTopic.loom import export_region_accessibility_to_loom, export_gene_activity_to_loom

def create_cluster_markers(diff_dict):
    return {
        'main': {
            cluster_name: df.index.tolist()
            for cluster_name, df in diff_dict.items()
        }
    }

# creat cluster markers
cluster_region_markers = create_cluster_markers(markers_dict)
cluster_gene_markers = create_cluster_markers(DAG_markers_dict)


os.makedirs(os.path.join(out_dir, "loom"), exist_ok=True)
export_region_accessibility_to_loom(
    accessibility_matrix = imputed_acc_obj,
    cistopic_obj = cistopic_obj,
    binarized_topic_region = region_bin_topics_otsu,
    binarized_cell_topic = binarized_cell_topic,
    selected_cells = cistopic_obj.projections['cell']['UMAP'].index.tolist(),
    out_fname = os.path.join(out_dir, "loom", "multiome_pycisTopic_region_accessibility.loom"),
    cluster_annotation = ['main'],
    cluster_markers = cluster_region_markers,
    tree_structure = ('multiome', 'pycisTopic'),
    title = 'multiome - Region accessibility all',
    nomenclature = "cattle",
    split_pattern = '-'
)
export_gene_activity_to_loom(
    gene_activity_matrix = 3,
    cistopic_obj = cistopic_obj,
    out_fname = os.path.join(out_dir, "loom", "multiome_pycisTopic_gene_activity.loom"),
    cluster_annotation = ['main'],
    cluster_markers = cluster_gene_markers,
    tree_structure = ('multiome', 'pycisTopic', 'ATAC'),
    title = 'multiome - Gene activity',
    nomenclature = "mm10",
    split_pattern = '-'
)
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_final.pkl"), "wb")
)
