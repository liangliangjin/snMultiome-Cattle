import pandas as pd
import numpy as np
import os

celltype_list = ['Nerve-cell-Neuron', 'Endothelial-cell-Cardiac-endothelial-cell-1', 'Nerve-cell-Granule-cell', 'Nerve-cell-GABAergic-neuron', 'Nerve-cell-Excitatory-neuron', 'Immune-cell-Microglia', 'Nerve-cell-Astrocyte', 'Nerve-cell-Oligodendrocyte-precursor-cell', 'Nerve-cell-Oligodendrocyte', 'Nerve-cell-Neural-precursor-cell', 'Stromal-cell-Stellate-cell', 'Unknown-cell', 'Epithelial-cell-Hepatocyte', 'Muscle-cell-Other-smooth-muscle-cell', 'Stromal-cell-Ruminal-fibroblast', 'Stromal-cell-Mesenchymal-cell', 'Epithelial-cell-Follicular-cell', 'Endothelial-cell-Ruminal-endothelial-cell-1', 'Immune-cell-Splenic-macrophage', 'Immune-cell-Cardiac-macrophage', 'Muscle-cell-Cardiomyocyte-1', 'Muscle-cell-Cardiomyocyte-2', 'Epithelial-cell-Alveolar-epithelial-type-II-cell', 'Immune-cell-NK-T-cell', 'Stromal-cell-Adipocyte-precursor-cell', 'Stromal-cell-Cardiac-fibroblast', 'Endothelial-cell-Liver-sinusoidal-endothelial-cell', 'Stromal-cell-Mesangial-cell', 'Endothelial-cell-Ovarian-endothelial-cell', 'Nerve-cell-Schwann-cell', 'Immune-cell-T-cell', 'Muscle-cell-Type-II-myonuclei', 'Endothelial-cell-Other-endothelial-cell', 'Immune-cell-Macrophage', 'Endothelial-cell-Cardiac-endothelial-cell-2', 'Stromal-cell-Satellite-cell', 'Muscle-cell-Ruminal-smooth-muscle-cell', 'Endothelial-cell-Splenic-sinusoidal-endothelial-cell', 'Stromal-cell-Splenic-fibroblast', 'Stromal-cell-Adipocyte', 'Muscle-cell-Cardiac-smooth-muscle-cell', 'Epithelial-cell-Keratinocyte', 'Epithelial-cell-Supporting-cell', 'Epithelial-cell-Proximal-convoluted-tubule-2', 'Muscle-cell-Type-I-myonuclei', 'Epithelial-cell-Proximal-convoluted-tubule-1', 'Epithelial-cell-Renal-tubular-epithelial-cell', 'Epithelial-cell-Loop-of-Henle---Distal-convoluted-tubule', 'Muscle-cell-Perivascular-smooth-muscle-like-cell', 'Immune-cell-B-cell', 'Stromal-cell-Lung-pericyte', 'Epithelial-cell-Alveolar-epithelial-type-I-cell', 'Germline-cell-Spermatid', 'Immune-cell-Muscle-macrophage', 'Endothelial-cell-Myoendothelial-cell', 'Stromal-cell-Myotendinous-junction', 'Muscle-cell-Hybrid-skeletal-muscle-fibers', 'Epithelial-cell-Hair-follicle-cell', 'Stromal-cell-Dermal-fibroblast', 'Endothelial-cell-Adipose-endothelial-cell', 'Epithelial-cell-Ruminal-epithelial-cell', 'Endothelial-cell-Ruminal-endothelial-cell-2', 'Immune-cell-Kupffer-cell', 'Germline-cell-Spermatocyte', 'Germline-cell-Spermatogonia', 'Immune-cell-Testicular-T-cell']

leiqiong_Sample = ['brain', 'heart', 'kidney', 'lung', 'muscle','spleen0', 'skin0', 'rumen', 'fat0', 'liver0', 'ovary', 'testis']
mongolian_Sample = ['brain2', 'heart2', 'kidney2', 'lung2', 'muscle2', 'spleen2', 'skin2', 'rumen2', 'fat2', 'liver2']


for temp_celltype in celltype_list:
    Used_df = pd.DataFrame()
    print(f"Processing cell type: {temp_celltype}")
    for temp_sample in leiqiong_Sample:
        file_name = temp_sample + '_' + temp_celltype + '.txt'
        if os.path.exists(file_name):
            temp_df = pd.read_csv(file_name, sep='\t')
            Used_df[temp_sample + '_' + temp_celltype] = temp_df['x']
    if not Used_df.empty:
        output_file = temp_celltype + '_leiqiong.txt'
        Used_df.to_csv(output_file, sep='\t')
        print(f"Saved file: {output_file}")
    else:
        print(f"No data for {temp_celltype}, skipping save.")



for temp_celltype in celltype_list:
    Used_df = pd.DataFrame()
    print(f"Processing cell type: {temp_celltype}")
    for temp_sample in mongolian_Sample:
        file_name = temp_sample + '_' + temp_celltype + '.txt'
        if os.path.exists(file_name):
            temp_df = pd.read_csv(file_name, sep='\t')
            Used_df[temp_sample + '_' + temp_celltype] = temp_df['x']
    if not Used_df.empty:
        output_file = temp_celltype + '_mongolian.txt'
        Used_df.to_csv(output_file, sep='\t')
        print(f"Saved file: {output_file}")
    else:
        print(f"No data for {temp_celltype}, skipping save.")




leiqiong_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_leiqiong.txt'):
        temp_leiqiong = pd.read_csv(temp_celltype+'_leiqiong.txt', sep='\t', index_col=0)
        leiqiong_df = pd.concat([leiqiong_df, temp_leiqiong], axis=1)
        print(temp_celltype)

leiqiong_df.to_csv('Settings/leiqiong_df.txt', sep='\t')


mongolian_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_mongolian.txt'):
        temp_mongolian = pd.read_csv(temp_celltype+'_mongolian.txt', sep='\t', index_col=0)
        mongolian_df = pd.concat([mongolian_df, temp_mongolian], axis=1)
        print(temp_celltype)

mongolian_df.to_csv('Settings/mongolian_df.txt', sep='\t')


for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_leiqiong.txt'):
        if os.path.exists(temp_celltype+'_mongolian.txt'):
            temp_leiqiong = pd.read_csv(temp_celltype+'_leiqiong.txt', sep='\t', index_col=0)
            temp_mongolian = pd.read_csv(temp_celltype+'_mongolian.txt', sep='\t', index_col=0)
            temp_df = pd.concat([temp_leiqiong, temp_mongolian], axis=1)
            temp_df.to_csv('Settings/'+temp_celltype+'_Cross.txt', sep='\t')

