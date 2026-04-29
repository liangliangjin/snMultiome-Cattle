import pandas as pd
import numpy as np
import os
import re

celltype_list = pd.read_csv('celltype_list', header=None)[0].tolist()

##DACRE, each breed and celltype
path = os.path.expanduser('~/SC/Peak_df/DA_CRE')
os.chdir(path)

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



##
os.makedirs('./Settings', exist_ok=True)
leiqiong_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_leiqiong.txt'):
        temp_leiqiong = pd.read_csv(temp_celltype+'_leiqiong.txt', sep='\t', index_col=0)
        leiqiong_df = pd.concat([leiqiong_df, temp_leiqiong], axis=1)
        print(temp_celltype)

leiqiong_df.to_csv('./Settings/leiqiong_df.txt', sep='\t')


mongolian_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_mongolian.txt'):
        temp_mongolian = pd.read_csv(temp_celltype+'_mongolian.txt', sep='\t', index_col=0)
        mongolian_df = pd.concat([mongolian_df, temp_mongolian], axis=1)
        print(temp_celltype)

mongolian_df.to_csv('./Settings/mongolian_df.txt', sep='\t')


for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'_leiqiong.txt'):
        if os.path.exists(temp_celltype+'_mongolian.txt'):
            temp_leiqiong = pd.read_csv(temp_celltype+'_leiqiong.txt', sep='\t', index_col=0)
            temp_mongolian = pd.read_csv(temp_celltype+'_mongolian.txt', sep='\t', index_col=0)
            temp_df = pd.concat([temp_leiqiong, temp_mongolian], axis=1)
            temp_df.to_csv('./Settings/'+temp_celltype+'_Cross.txt', sep='\t')


##Only celltype
path = os.path.expanduser('~/SC/Peak_df/DA_CRE_Onlycelltypes')
os.chdir(path)
cattle_df = pd.DataFrame()
for temp_celltype in celltype_list:
    if os.path.exists(temp_celltype+'.txt'):
        temp_cattle = pd.read_csv(temp_celltype+'.txt', sep='\t', index_col=0)
        temp_cattle.columns = [temp_celltype] * temp_cattle.shape[1]
        cattle_df = pd.concat([cattle_df, temp_cattle], axis=1)
        print(temp_celltype)

cattle_df.to_csv('cattle_df.txt', sep='\t')

##tissue and celltype
path = os.path.expanduser('~/SC/Peak_df/DA_CRE_byTissue')
os.chdir(path)
output_dir = './Settings'
os.makedirs(output_dir, exist_ok=True)

all_df = pd.DataFrame()
meta_list = []

all_files = [f for f in os.listdir(path) if f.endswith('.txt')]

for file_name in all_files:
    file_path = os.path.join(path, file_name)
    base = file_name[:-4]
    matched = False
    for ct in celltype_list:
        suffix = '_' + ct
        if base.endswith(suffix):
            tissue = base[:-len(suffix)]
            celltype = ct
            matched = True
            break
    if not matched:
        print(f'Skip: cannot parse {file_name}')
        continue
    temp_df = pd.read_csv(file_path, sep='\t', index_col=0)
    temp_df.columns = [base]
    all_df = pd.concat([all_df, temp_df], axis=1)
    meta_list.append({
        'sample_id': base,
        'tissue': tissue,
        'celltype': celltype
    })
    print(f'Loaded: {base}')

all_df.to_csv(os.path.join(output_dir, 'tissue_celltype_matrix.txt'), sep='\t')

meta_df = pd.DataFrame(meta_list)
meta_df.to_csv(os.path.join(output_dir, 'tissue_celltype_meta.txt'), sep='\t', index=False)

for temp_celltype in celltype_list:
    ct_cols = meta_df.loc[meta_df['celltype'] == temp_celltype, 'sample_id'].tolist()
    if len(ct_cols) > 0:
        temp_cross = all_df[ct_cols]
        temp_cross.to_csv(os.path.join(output_dir, f'{temp_celltype}_Cross.txt'), sep='\t')
        print(f'Saved: {temp_celltype}_Cross.txt')