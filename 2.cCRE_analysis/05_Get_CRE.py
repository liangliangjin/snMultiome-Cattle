import pandas as pd
import numpy as np
import os

celltype_list = ['Nerve-cell-Neuron', 'Endothelial-cell-Cardiac-endothelial-cell-1', 'Nerve-cell-Granule-cell', 'Nerve-cell-GABAergic-neuron', 'Nerve-cell-Excitatory-neuron', 'Immune-cell-Microglia', 'Nerve-cell-Astrocyte', 'Nerve-cell-Oligodendrocyte-precursor-cell', 'Nerve-cell-Oligodendrocyte', 'Nerve-cell-Neural-precursor-cell', 'Stromal-cell-Stellate-cell', 'Unknown-cell', 'Epithelial-cell-Hepatocyte', 'Muscle-cell-Other-smooth-muscle-cell', 'Stromal-cell-Ruminal-fibroblast', 'Stromal-cell-Mesenchymal-cell', 'Epithelial-cell-Follicular-cell', 'Endothelial-cell-Ruminal-endothelial-cell-1', 'Immune-cell-Splenic-macrophage', 'Immune-cell-Cardiac-macrophage', 'Muscle-cell-Cardiomyocyte-1', 'Muscle-cell-Cardiomyocyte-2', 'Epithelial-cell-Alveolar-epithelial-type-II-cell', 'Immune-cell-NK-T-cell', 'Stromal-cell-Adipocyte-precursor-cell', 'Stromal-cell-Cardiac-fibroblast', 'Endothelial-cell-Liver-sinusoidal-endothelial-cell', 'Stromal-cell-Mesangial-cell', 'Endothelial-cell-Ovarian-endothelial-cell', 'Nerve-cell-Schwann-cell', 'Immune-cell-T-cell', 'Muscle-cell-Type-II-myonuclei', 'Endothelial-cell-Other-endothelial-cell', 'Immune-cell-Macrophage', 'Endothelial-cell-Cardiac-endothelial-cell-2', 'Stromal-cell-Satellite-cell', 'Muscle-cell-Ruminal-smooth-muscle-cell', 'Endothelial-cell-Splenic-sinusoidal-endothelial-cell', 'Stromal-cell-Splenic-fibroblast', 'Stromal-cell-Adipocyte', 'Muscle-cell-Cardiac-smooth-muscle-cell', 'Epithelial-cell-Keratinocyte', 'Epithelial-cell-Supporting-cell', 'Epithelial-cell-Proximal-convoluted-tubule-2', 'Muscle-cell-Type-I-myonuclei', 'Epithelial-cell-Proximal-convoluted-tubule-1', 'Epithelial-cell-Renal-tubular-epithelial-cell', 'Epithelial-cell-Loop-of-Henle---Distal-convoluted-tubule', 'Muscle-cell-Perivascular-smooth-muscle-like-cell', 'Immune-cell-B-cell', 'Stromal-cell-Lung-pericyte', 'Epithelial-cell-Alveolar-epithelial-type-I-cell', 'Germline-cell-Spermatid', 'Immune-cell-Muscle-macrophage', 'Endothelial-cell-Myoendothelial-cell', 'Stromal-cell-Myotendinous-junction', 'Muscle-cell-Hybrid-skeletal-muscle-fibers', 'Epithelial-cell-Hair-follicle-cell', 'Stromal-cell-Dermal-fibroblast', 'Endothelial-cell-Adipose-endothelial-cell', 'Epithelial-cell-Ruminal-epithelial-cell', 'Endothelial-cell-Ruminal-endothelial-cell-2', 'Immune-cell-Kupffer-cell', 'Germline-cell-Spermatocyte', 'Germline-cell-Spermatogonia', 'Immune-cell-Testicular-T-cell']

# Get cell type specific cCRE Set (leiqiong)
for temp_celltype in celltype_list:
    if os.path.exists('./Res/'+temp_celltype+'_leiqiong.txt'):
        cattle_df = pd.read_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE/Res/'+temp_celltype+'_leiqiong.txt', sep='\t')
        cattle_df = cattle_df.loc[cattle_df['padj']<0.05,]
        cattle_df = cattle_df.loc[cattle_df['logFC']<0,]
        if os.path.exists('/home/Jingliangliang/SC/Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt'):
            Cross_df = pd.read_csv('Res/'+temp_celltype+'_Cross_IvsT.txt', sep='\t')
            Cross_df = Cross_df.loc[Cross_df['padj']<0.05,]
            Cross_df = Cross_df.loc[Cross_df['logFC']<0,]
            print(temp_celltype, cattle_df.shape, Cross_df.shape)
            Used_peak = np.intersect1d(cattle_df.index, Cross_df.index)
            Used_df = cattle_df.loc[Used_peak, ]
            Used_df = Used_df.sort_values(by='padj')
            print(Used_df.shape)
            Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
            Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
            Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]
            if len(Used_df)>0:
                Used_df.loc[:, ['chr','start','end','logFC','padj']].to_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE/DAcCRE_celltype/'+temp_celltype+'_leiqiong.txt', sep='\t', index=None, header=None)


# Get cell type specific cCRE Set (mongolian)
for temp_celltype in celltype_list:
    if os.path.exists('./Res/'+temp_celltype+'_mongolian.txt'):
        cattle_df = pd.read_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE/Res/'+temp_celltype+'_mongolian.txt', sep='\t')
        cattle_df = cattle_df.loc[cattle_df['padj']<0.05,]
        cattle_df = cattle_df.loc[cattle_df['logFC']<0,]
        if os.path.exists('/home/Jingliangliang/SC/Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt'):
            Cross_df = pd.read_csv('Res/'+temp_celltype+'_Cross_IvsT.txt', sep='\t')
            Cross_df = Cross_df.loc[Cross_df['padj']<0.05,]
            Cross_df = Cross_df.loc[Cross_df['logFC']>0,]
            print(temp_celltype, cattle_df.shape, Cross_df.shape)
            Used_peak = np.intersect1d(cattle_df.index, Cross_df.index)
            Used_df = cattle_df.loc[Used_peak, ]
            Used_df = Used_df.sort_values(by='padj')
            print(Used_df.shape)
            Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
            Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
            Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]
            if len(Used_df)>0:
                Used_df.loc[:, ['chr','start','end','logFC','padj']].to_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE/DAcCRE_celltype/'+temp_celltype+'_mongolian.txt', sep='\t', index=None, header=None)



# Get cell type specific cCRE Set (cattle)
for temp_celltype in celltype_list:
    if os.path.exists('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/Res/'+temp_celltype+'.txt'):
        cattle_df = pd.read_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/Res/'+temp_celltype+'.txt', sep='\t')
        cattle_df = cattle_df.loc[cattle_df['padj']<0.05,]
        cattle_df = cattle_df.loc[cattle_df['logFC']<0,]
        print(temp_celltype, cattle_df.shape)
        Used_df = cattle_df
        Used_df = Used_df.sort_values(by='padj')
        print(Used_df.shape)
        Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
        Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
        Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]
        if len(Used_df)>0:
            Used_df.loc[:, ['chr','start','end','logFC','padj']].to_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/Res/P005/'+temp_celltype+'_P005.txt', sep='\t', index=None, header=None)

final_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/Res/'+temp_celltype+'.txt'):
        cattle_df = pd.read_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/Res/'+temp_celltype+'.txt', sep='\t')
        cattle_df['value'] = ((cattle_df['padj'] < 0.05) & (cattle_df['logFC'] < 0)).astype(int)
        column_name = temp_celltype
        final_df[column_name] = cattle_df['value']
        count_1 = (cattle_df['value'] == 1).sum()
        print(f'{temp_celltype}: {count_1} rows where padj < 0.05 and logFC < 0')

final_df['chr'] = [x.split('_')[0] for x in cattle_df.index]
final_df['start'] = [int(x.split('_')[1]) for x in cattle_df.index]
final_df['end'] = [int(x.split('_')[2]) for x in cattle_df.index]
final_df = final_df[['chr', 'start', 'end'] + [col for col in final_df.columns if col not in ['chr', 'start', 'end']]]
final_df.to_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE_Onlycelltypes/merged_results.txt', sep='\t', index=False)


