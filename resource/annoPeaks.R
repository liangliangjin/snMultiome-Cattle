#Referring to the annotation function of ArchR, more peaktypes have been added.
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(org.Bt.eg.db)
library(ArchR)
addArchRThreads(threads = 16)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
txdb <- makeTxDbFromGFF("ARS-UCD1.2.110.chr.gtf", format = "gtf") # Or replace it with the GTF file of another species/version.
geneAnnotation <- createGeneAnnotation(TxDb = txdb, OrgDb = org.Bt.eg.db)
geneAnnotation$cds <- unlist(cdsBy(txdb, "tx"))
geneAnnotation$fiveUTRs <- unlist(fiveUTRsByTranscript(txdb))
geneAnnotation$threeUTRs <- unlist(threeUTRsByTranscript(txdb))
geneAnnotation$cds$gene_id <- mapIds(txdb, names(geneAnnotation$cds), "GENEID", "TXID")
geneAnnotation$fiveUTRs$gene_id <- mapIds(txdb, names(geneAnnotation$fiveUTRs), "GENEID", "TXID")
geneAnnotation$threeUTRs$gene_id <- mapIds(txdb, names(geneAnnotation$threeUTRs), "GENEID", "TXID")
geneAnnotation$cds$symbol <- geneAnnotation$genes$symbol[match(geneAnnotation$cds$gene_id, geneAnnotation$genes$gene_id)]
geneAnnotation$fiveUTRs$symbol <- geneAnnotation$genes$symbol[match(geneAnnotation$fiveUTRs$gene_id, geneAnnotation$genes$gene_id)]
geneAnnotation$threeUTRs$symbol <- geneAnnotation$genes$symbol[match(geneAnnotation$threeUTRs$gene_id, geneAnnotation$genes$gene_id)]
names(geneAnnotation)
if(!any(grepl("^chr", seqlevels(geneAnnotation$genes)))){
	seqlevels(geneAnnotation$genes) <- paste0("chr", seqlevels(geneAnnotation$genes))
	seqlevels(geneAnnotation$exons) <- paste0("chr", seqlevels(geneAnnotation$exons))
	seqlevels(geneAnnotation$TSS) <- paste0("chr", seqlevels(geneAnnotation$TSS))
	seqlevels(geneAnnotation$cds) <- paste0("chr", seqlevels(geneAnnotation$cds))
	seqlevels(geneAnnotation$fiveUTRs) <- paste0("chr", seqlevels(geneAnnotation$fiveUTRs))
	seqlevels(geneAnnotation$threeUTRs) <- paste0("chr", seqlevels(geneAnnotation$threeUTRs))
}
validGRanges <- function(gr = NULL){
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}
BSgenome = eval(parse(text = genomeAnnotation$genome))
annoPeaks <- function(peaks){
	peaks <- validGRanges(peaks)
	peakSummits <- GenomicRanges::resize(peaks,1,"center")
	geneAnnotation$genes <- geneAnnotation$genes
	geneAnnotation$exons <- geneAnnotation$exons
	geneAnnotation$TSS <- geneAnnotation$TSS
	geneAnnotation$fiveUTRs <- geneAnnotation$fiveUTRs
    geneAnnotation$threeUTRs <- geneAnnotation$threeUTRs
	BSgenome <- validBSgenome(BSgenome)
	distPeaks <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
	mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
	# The defined region of the promoter can be modified to be consistent with other downstream analyses.
	promoters <- extendGR(GenomicRanges::resize(geneAnnotation$genes, 1, "start"), upstream = 1000, downstream = 1000)
	op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
	og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
	oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
	oc <- overlapsAny(peakSummits, geneAnnotation$cds, ignore.strand = TRUE)
	of <- overlapsAny(peakSummits, geneAnnotation$fiveUTRs, ignore.strand = TRUE)
	ot <- overlapsAny(peakSummits, geneAnnotation$threeUTRs, ignore.strand = TRUE)
	peakType_simple <- rep("Distal", length(peaks))
	peakType_simple[which(og & oe)] <- "Exonic"
	peakType_simple[which(og & !oe)] <- "Intronic"
	peakType_simple[which(op)] <- "Promoter"
	mcols(peaks)$peakType <- peakType_simple
	
	type_matrix <- cbind(
		`5'UTR` = of,
		`3'UTR` = ot,
		CDS = oc,
		Intronic = og & !oe,
		Exonic = oe & !oc & !of & !ot,
		Promoter = op,
		Distal = !og & !op
	)
	type_list <- apply(type_matrix, 1, function(x) {
		types <- colnames(type_matrix)[x]
		if(length(types) == 0) return("Distal")
		paste(types, collapse = ",")
	})
	# If the peak covers multiple type areas, display all of them and separate them with commas.
	mcols(peaks)$peakTypeMulti <- type_list
	distTSS <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToTSS <- mcols(distTSS)$distance
	if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distTSS)]
	}else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distTSS)]
	}
	nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  	mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  	mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  	return(peaks)
}
## The input file should be in Grange format.
#peaks_gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
# Usage: peaks_anno <- annoPeaks(peaks_gr).