###########
#EdgeR process - only celltypes is considered
###########
library(edgeR)
library(stringr)

output_dir <- './Peak_df/DA_CRE_Onlycelltypes/Res/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

celltype_list <- readLines("celltype_list")
Used_order = c()

for (temp_celltype in celltype_list) {
    if (file.exists(paste0('./Peak_df/DA_CRE_Onlycelltypes/',temp_celltype, '.txt'))) {
        Used_order <- c(Used_order, temp_celltype)
    }
}
bcv = 0.1
Count_df = read.delim('./Peak_df/DA_CRE_Onlycelltypes/cattle_df.txt', sep='\t', check.names = FALSE, row.names = 1)
for (temp_celltype in celltype_list){
    Used_celltype = c()
    for (it in Used_order){
        if (it == temp_celltype){
            Used_celltype = c(Used_celltype, temp_celltype)
        }
        else{
            Used_celltype = c(Used_celltype, 'Other')
        }
    }
    Count_df_sub = Count_df[, colnames(Count_df) %in% Used_order]
    meta = data.frame(Celltype=Used_celltype, row.names = Used_order)
    if(length(unique(meta$Celltype)) > 1){
		message("Test: ", temp_celltype)
		data.use = as.matrix(Count_df_sub)
		group <- factor(meta$Celltype, levels = c(temp_celltype, "Other"))
		y <- DGEList(counts = data.use, group = group)
		y <- calcNormFactors(y)
		tb.pos <- exactTest(y, dispersion = bcv^2)$table;
		tb.pos$padj = p.adjust(tb.pos$PValue, method = 'BH')
		output_file = paste0(output_dir, temp_celltype, '.txt')
		write.table(tb.pos, output_file, sep='\t', quote = FALSE)
		message("save ", output_file)
	}
}



###########
#EdgeR process for Leiqiong and mongolian cattle
###########
library(edgeR)
library(stringr)
celltype_list <- readLines("celltype_list")

leiqiong_Sample = c('brain', 'heart', 'kidney', 'lung', 'muscle','spleen0', 'skin0', 'rumen', 'fat0', 'liver0', 'ovary', 'testis')
mongolian_Sample = c('brain2', 'heart2', 'kidney2', 'lung2', 'muscle2','spleen2', 'skin2', 'rumen2', 'fat2', 'liver2')
bcv = 0.1
dir.create("./Peak_df/DA_CRE/Res", showWarnings = FALSE)
for (breed in c("leiqiong", "mongolian")) {
    if (breed == "leiqiong") {
        sample_list <- leiqiong_Sample
        count_file <- './Peak_df/DA_CRE/Settings/leiqiong_df.txt'
        output_prefix <- '_leiqiong.txt'
    } else {
        sample_list <- mongolian_Sample
        count_file <- './Peak_df/DA_CRE/Settings/mongolian_df.txt'
        output_prefix <- '_mongolian.txt'
    }
    Used_order = c()
    for (temp_celltype in celltype_list) {
        for (temp_sample in sample_list) {
            if (file.exists(paste0('./Peak_df/DA_CRE/', temp_sample, '_', temp_celltype, '.txt'))) {
                Used_order <- c(Used_order, paste0(temp_sample, '_', temp_celltype))
            }
        }
    }
    Count_df = read.delim(count_file, sep='\t', check.names = FALSE, row.names=1)
    for (temp_celltype in celltype_list){
        Used_celltype = c()
        for (it in Used_order){
            parts <- strsplit(it, '_')[[1]]
            current_celltype <- parts[length(parts)]
            if (current_celltype == temp_celltype){
                Used_celltype = c(Used_celltype, temp_celltype)
            } else {
                Used_celltype = c(Used_celltype, 'Other')
            }
        }
        Count_df_sub <- Count_df[, colnames(Count_df) %in% Used_order, drop = FALSE]
        meta = data.frame(Celltype = Used_celltype, row.names = Used_order)
        if(length(unique(meta$Celltype)) > 1){
            message("Testing: ", temp_celltype, " in ", breed)
            data.use = as.matrix(Count_df_sub)
            group <- factor(meta$Celltype, levels = c(temp_celltype, "Other"))
            y <- DGEList(counts = data.use, group = group)
            y <- calcNormFactors(y)
            tb.pos <- exactTest(y, dispersion = bcv^2)$table
            tb.pos$padj = p.adjust(tb.pos$PValue, method = 'BH')
            output_file = paste0('./Peak_df/DA_CRE/Res/', temp_celltype, output_prefix)
            write.table(tb.pos, output_file, sep='\t', quote = FALSE)
            message("Saved ", output_file)
        }
    }
}

# edgeR test for cross-subspecies comparision
bcv = 0.1
for (temp_celltype in celltype_list){
    Used_order = c()
    for (it in leiqiong_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    for (it in mongolian_Sample){
        Used_order = c(Used_order, paste0(it, '_', temp_celltype))
    }
    Used_celltype = c(rep('leiqiong', length(leiqiong_Sample)), 
                      rep('mongolian', length(mongolian_Sample)))
    meta = data.frame(Celltype=Used_celltype, row.names=Used_order)
    input_file = paste0('./Peak_df/DA_CRE/Settings/', temp_celltype, '_Cross.txt')
    output_file = paste0('./Peak_df/DA_CRE/Res/', temp_celltype, '_Cross_IvsT.txt')
    if(file.exists(input_file)){
		Count_df = read.delim(input_file, sep='\t', check.names = FALSE, row.names=1)
		Count_df = Count_df[, gsub("\\.", "-", colnames(Count_df))%in%Used_order]
		data.use = as.matrix(Count_df)
		group <- factor(meta$Celltype[Used_order%in%colnames(Count_df)], levels = c("leiqiong", "mongolian"))
		y <- DGEList(counts=data.use, group=group)
		y <- calcNormFactors(y)
		tb.pos <- exactTest(y, dispersion=bcv^2)$table;
		tb.pos$padj = p.adjust(tb.pos$PValue, method = 'BH')
		write.table(tb.pos, output_file, sep='\t', quote=FALSE)
		message(temp_celltype)
    }
}