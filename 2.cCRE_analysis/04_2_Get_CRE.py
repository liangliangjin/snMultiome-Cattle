import pandas as pd
import numpy as np
import os

celltype_list = pd.read_csv('celltype_list', header=None)[0].tolist()
os.makedirs('./Peak_df/DA_CRE/DAcCRE_celltype', exist_ok=True)
# Get cell type specific cCRE Set (leiqiong)
for temp_celltype in celltype_list:
    if os.path.exists('./Peak_df/DA_CRE/Res/'+temp_celltype+'_leiqiong.txt'):
        cattle_df = pd.read_csv('./Peak_df/DA_CRE/Res/'+temp_celltype+'_leiqiong.txt', sep='\t')
        cattle_df = cattle_df.loc[cattle_df['padj'] < 0.01,]
        cattle_df = cattle_df.loc[cattle_df['logFC']< 0,]
        if os.path.exists('./Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt'):
            Cross_df = pd.read_csv('./Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt', sep='\t')
            Cross_df = Cross_df.loc[Cross_df['padj']< 0.01,]
            Cross_df = Cross_df.loc[Cross_df['logFC']< 0,]
            print(temp_celltype, cattle_df.shape, Cross_df.shape)
            Used_peak = np.intersect1d(cattle_df.index, Cross_df.index)
            Used_df = cattle_df.loc[Used_peak, ]
            Used_df = Used_df.sort_values(by='padj')
            print(Used_df.shape)
            Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
            Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
            Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]
            if len(Used_df)>0:
                Used_df.loc[:, ['chr','start','end','logFC','padj']].to_csv('./Peak_df/DA_CRE/DAcCRE_celltype/'+temp_celltype+'_leiqiong.txt', sep='\t', index=None, header=None)


# Get cell type specific cCRE Set (mongolian)
for temp_celltype in celltype_list:
    if os.path.exists('./Peak_df/DA_CRE/Res/'+temp_celltype+'_mongolian.txt'):
        cattle_df = pd.read_csv('/home/Jingliangliang/SC/Peak_df/DA_CRE/Res/'+temp_celltype+'_mongolian.txt', sep='\t')
        cattle_df = cattle_df.loc[cattle_df['padj'] < 0.01,]
        cattle_df = cattle_df.loc[cattle_df['logFC'] < 0,]
        if os.path.exists('./Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt'):
            Cross_df = pd.read_csv('./Peak_df/DA_CRE/Res/'+temp_celltype+'_Cross_IvsT.txt', sep='\t')
            Cross_df = Cross_df.loc[Cross_df['padj'] < 0.01,]
            Cross_df = Cross_df.loc[Cross_df['logFC'] > 0,]
            print(temp_celltype, cattle_df.shape, Cross_df.shape)
            Used_peak = np.intersect1d(cattle_df.index, Cross_df.index)
            Used_df = cattle_df.loc[Used_peak, ]
            Used_df = Used_df.sort_values(by='padj')
            print(Used_df.shape)
            Used_df['chr'] = [x.split('_')[0] for x in Used_df.index]
            Used_df['start'] = [int(x.split('_')[1]) for x in Used_df.index]
            Used_df['end'] = [int(x.split('_')[2]) for x in Used_df.index]
            if len(Used_df)>0:
                Used_df.loc[:, ['chr','start','end','logFC','padj']].to_csv('./Peak_df/DA_CRE/DAcCRE_celltype/'+temp_celltype+'_mongolian.txt', sep='\t', index=None, header=None)



# #merge.sh
#!/bin/bash                                                                                                                                                                         
# if [ $# -ne 1 ]; then
	# echo "Usage: $0 Directory path"
	# exit 1
# fi

# dir=$1
# output="CRE_merge.txt"

# rm -f "$output"

# for file in "$dir"/*.txt; do
	# filename=$(basename "$file" .txt)
	# awk -v fname="$filename" '{print $0 "\t" fname}' "$file" >> "$output"
# done
# echo "The merge is complete, and the result is saved in $output"
