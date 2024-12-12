library(edgeR)
library(stringr)
celltype_list <- readLines("celltype_list")
Used_order = c()
leiqiong_Sample = c('brain', 'heart', 'kidney', 'lung', 'muscle','spleen0', 'skin0', 'rumen', 'fat0', 'liver0', 'ovary', 'testis')
for (temp_celltype in celltype_list) {
    for (temp_sample in leiqiong_Sample) {
        if (file.exists(paste0('/DA_CRE/',temp_sample, '_', temp_celltype, '.txt'))) {
            Used_order <- c(Used_order, paste0(temp_sample, '_', temp_celltype))
        }
    }
}
bcv = 0.1
for (temp_celltype in celltype_list){
    Count_df = read.delim('/DA_CRE/Settings/leiqiong_df.txt', sep='\t', row.names=1)
    Used_celltype = c()
    for (it in Used_order){
        if (strsplit(it, '_')[[1]][2]==temp_celltype){
            Used_celltype = c(Used_celltype, temp_celltype)
        }
        else{
            Used_celltype = c(Used_celltype, 'Other')
        }
    }
    Count_df = Count_df[, gsub("\\.","-",colnames(Count_df))%in%Used_order]
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)
    if(length(unique(meta$Celltype))>1){
    data.use = as.matrix(Count_df)
    group <- factor(meta$Celltype)
    design <- model.matrix(~group)
    y <- DGEList(counts=data.use, group=group)
    tb.pos <- exactTest(y, dispersion=bcv^2)$table;
    tb.pos$padj = p.adjust(tb.pos$PValue, method = 'hochberg')
    output_file = paste0('/DA_CRE/Res/', temp_celltype, '_leiqiong.txt')
    write.table(tb.pos, output_file, sep='\t', quote=FALSE)
	print(paste0("save ",output_file))
	}
}
###########
#Run the same edgeR process for Mongolian cattle
###########

# edgeR test for cross-species comparision
bcv = 0.1
for (temp_celltype in celltype_list){
    Used_order = c()
    for (it in leiqiong_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    for (it in mongolian_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    Used_celltype = c(rep('leiqiong',length(leiqiong_Sample)), 
                      rep('mongolian',length(mongolian_Sample)))
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)
    input_file = paste0('/DA_CRE/Settings/', temp_celltype, '_Cross.txt')
    output_file = paste0('/DA_CRE/Res/', temp_celltype, '_Cross_IvsT.txt')
    if(file.exists(input_file)){
    Count_df = read.delim(input_file, sep='\t', row.names=1)
    Count_df = Count_df[, gsub("\\.","-",colnames(Count_df))%in%Used_order]
    data.use = as.matrix(Count_df)
    group <- factor(meta$Celltype[Used_order%in%gsub("\\.","-",colnames(Count_df))])
    design <- model.matrix(~group)
    y <- DGEList(counts=data.use, group=group)
    tb.pos <- exactTest(y, dispersion=bcv^2)$table;
    tb.pos$padj = p.adjust(tb.pos$PValue, method = 'hochberg')
    write.table(tb.pos, output_file, sep='\t', quote=FALSE)
	message(temp_celltype)
    }
}