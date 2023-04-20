library(rtracklayer)

setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code")

# Data from Zhang2021 A single cell atlas of chromatin accessibility in the human genome
# * 30 adult and 15 fetal tissues, scATAC-seq
# * 1,152,611 cCREs across 111(adult)+111(fetal)=222 cell types
#   * 890,130 cCREs in adult cell types
#    * 802,025 cCREs in fetal cell types

# cCRE_hg38.tsv.gz
# ================
#
#   The list of 1,154,611 cCREs. There are 7 columns in this file:
#
# 1. chromosome: name of the chromosome.
# 2. hg38_Start: 0-based starting location of the cCRE in the genome (hg38).
# 3. hg38_End: End of the cCRE in the genome (hg38).
# 4. class: Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
# 5. Present in fetal tissues: if this cCRE is detected in at least one fetal tissue
# 6. Present in adult tissues: if this cCRE is detected in at least one adult tissue
# 7. CRE module: The ID of CRE module that the cCRE belongs to. 1-150

#What is the definition of cCRE and CRE module
cCRE_hg38=read.table("../CATLAS/cCRE_hg38.tsv.gz",sep="\t", header=FALSE)

Cell_ontology=read.table("../CATLAS/Cell_ontology.tsv",sep="\t", header=TRUE)

#Cell_metadata.tsv.gz
# ====================
#   
#   This file contains metadata for each cell/barcode. It has 6 columns:
#   
# 1. cellID: unique ID.
# 2. logUMI: log10 number of unique fragment passing QC.
# 3. tsse: TSS enrichment
# 4. tissue: The ID of tissue sample
# 5. cell type: The name of the cell type
# 6. Life stage: Fetal or Adult

Cell_metadata=read.table("../CATLAS/Cell_metadata.tsv.gz",sep="\t", header=TRUE)

Celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t", header=FALSE)

# * 435,142 adult cCREs with restricted accessibility in one or a few cell types


#From Mendeley database
files=dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/")

ct_names=do.call(rbind,sapply(files, strsplit, ".bed.gz"))
dimnames(ct_names)=NULL

cCREs_list<-GRangesList()

for(ct in ct_names){
  print(ct)
  cCREs_list[[ct]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/",ct,".bed.gz"),format="bed")


}

saveRDS(cCREs_list,  file = paste0( "../RData/cCREs_list.Rds")) #
#cCREs_list <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/cCREs_list.Rds")

sum(sapply(cCREs_list, length)) #1,091,805

#unique
length(unique(unlist(cCREs_list))) #435142 adults cCREs with restricted accessibility in one or a few cell types


# * Phylogenetic analysis to group fetal and adult cell types into different cell lineages, enrichment for lineage-specific TFBSs?
# * Differentially accessible (DA) elements between adults and fetal cell lines/cell lineage groups.
#  Adult and fetal-specific cCREs for the cell-lineage groups. Specific cCREs underlying fetal or adult-specific regulatory programs.

#WHERE IS THIS DATA?

summary=read.csv("/scratch/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/yv4fzv6cnm-4/Zhang et al Figure 4/Diff_peaks_motifs/summary.tsv", sep="\t", header=TRUE)

cell_lineages=summary$Group.name #17

dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 4/Diff_peaks_motifs/") #same as cell_lineages

adult_cCREs<-GRangesList()
fetal_cCREs<-GRangesList()

for(cl in cell_lineages){
  print(cl)
  adult_cCREs[[cl]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 4/Diff_peaks_motifs/",cl,"/adult.bed"),format="bed") #30154
  fetal_cCREs[[cl]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 4/Diff_peaks_motifs/",cl,"/fetal.bed"),format="bed") #25958
  
}

saveRDS(adult_cCREs,  file = paste0( "../RData/adult_cCREs.Rds")) #
saveRDS(fetal_cCREs,  file = paste0( "../RData/fetal_cCREs.Rds")) #


#cCRE_hg38=read.table("../CATLAS/cCRE_hg38.tsv.gz",sep="\t", header=FALSE)
#nrow(cCRE_hg38)

#Cell_ontology=read.table("../CATLAS/Cell_ontology.tsv",sep="\t", header=TRUE)
# 1. chromosome: name of the chromosome.
# 2. hg38_Start: 0-based starting location of the cCRE in the genome (hg38).
# 3. hg38_End: End of the cCRE in the genome (hg38).
# 4. class: Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
# 5. Present in fetal tissues: if this cCRE is detected in at least one fetal tissue
# 6. Present in adult tissues: if this cCRE is detected in at least one adult tissue
# 7. CRE module: The ID of CRE module that the cCRE belongs to.

#1,152,611 cCREs across 111(adult)+111(fetal)=222 cell types

# * Cell-type specificity of cCREs across the cell types.
# 1.2 million cCREs clustered into 150 cis-regulatory modules (CRMs).
# Putative master regulators of fetal and adult cell types.


library(readr)
gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")

cCREs=read_delim("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_cCRE_list_GRCh38.tsv.gz", delim="\t", col_names=TRUE)
nrow(cCREs) #1154611, this is the same table as cCRE_hg38

names(cCREs)=c("chr", "start", "end", "class", "fetal", "adult", "CREmoduleID")

library("tidyverse")
library("GenomicRanges")
library("dplyr")



by_CREdmodules <- cCREs %>% group_by(CREmoduleID)

by_CREdmodules %>% tally(sort=TRUE)

group_indices=by_CREdmodules %>% group_indices()

CREmodules <- GRangesList()

for(gk in 1:150){
  print(gk)
  
  data=as.data.frame(cCREs %>% filter( CREmoduleID==gk) %>% select(c("chr", "start", "end")))
  tmp_GRanges=GRanges(seqnames = Rle( data$chr, rep(1, nrow(data)) ),
                      ranges = IRanges(start=data$start, end = data$end) ,
                      strand = Rle(strand( rep("*", nrow(data))) , rep(1, nrow(data))  ))
  
  
  isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  
  seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  
  #try dropSeqlevels or ?keepStandardChromosomes from GenomicRanges library.
  
  
  tmp_GRanges=dropSeqlevels(tmp_GRanges, names(seqlengths(tmp_GRanges)[ which( (as.vector(names(seqlengths(tmp_GRanges))) %in% names(seqlengths)==FALSE) ) ]),
                            pruning.mode = "coarse")
  
  #tmp_GRanges=tmp_GRanges[which( (as.vector(seqnames(tmp_GRanges)) %in% names(seqlengths)==FALSE) )]
  
  #seqlengths(tmp_GRanges)= seqlengths(tmp_GRanges)[ which( (as.vector(names(seqlengths(tmp_GRanges))) %in% names(seqlengths)==FALSE) ) ]
  
  #CREmodules may contain some chromosomes not in gtf
  
  
  seqlengths(tmp_GRanges)=seqlengths[names( seqlengths(tmp_GRanges))]
  
  genome(tmp_GRanges) <- "GRCh38"
  
  CREmodules[[gk]]=tmp_GRanges
  
  
}


#saveRDS(CREmodules,  file ="../RData/CREmodules.Rds") 
CREmodules<-readRDS(file = "../RData/CREmodules.Rds")

#Modules connected to cell_types, groups of cell types

Meta_modules=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5D_CRE_Module_Motif_enrichment.csv",sep="\t")
Meta_modules=as.data.frame(t(Meta_modules[c(1,2), -1]))
rownames(Meta_modules)=NULL
names(Meta_modules)=c("ModuleGroup", "ModuleID")
Meta_modules[,"ModuleID"]=as.numeric(Meta_modules[,"ModuleID"])

saveRDS(Meta_modules,  file ="../RData/Meta_modules.Rds") 




#qn=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Quantile_normalized_log_RPKM.tsv-NN6XGL.gz",sep="\t", header=TRUE)

RAW_RPKM=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Raw_RPKM.tsv.gz",sep="\t", header=TRUE)

cell_types=colnames(RAW_RPKM)[-1] #222

#RAW_RPKM=read_delim("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Raw_RPKM.tsv.gz", delim="\t", col_names=TRUE)
#RAW_RPKM[1:5,]


#Aims:

#Motif enrichment analysis on cell-type restricted cCREs
  #345,142 cCREs that demonstrated restricted accessibility in one or a few (how many?) cell types
  #to reveal putative TFs of each cell type
  #cCREs restricted to a certain cell type, 10 most significant TFs selected of a single cell type
  #results is cell type x TF matrix
  # Color represents -log_10(P-val), depleted (blue), enriched (red)


#cCRE_by_cell_type
#=================
#This directory stores cCRE by cell type, containing three files:

#1. matrix.mtx.gz: matrix in matrix market format, each line represents a cCRE-cell type pair. For example, "2 22" means cCRE No. 2 is accessibile in cell type No. 22.
#2. celltypes.txt.gz: column names (cell types) of the matrix.
#3. cCREs.bed.gz: row names (cCREs) of the matrix.
library("Matrix")
matrix=readMM("../CATLAS/matrix.tsv.gz") #this is 1,154,611 x 222 matrix, contains TRUE or FALSE

tmp=do.call(rbind, apply(matrix,1,table)) #For each cCREs, in how many it is present
restricted_cCREs=which(tmp[,"TRUE"]==1) #cCREs that are present in a single cell type 347991 (paper: 345142 cCREs, maybe blacklist removed)

celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t")


#computed the odds ratio and the significance of enrichment in each cluster,
#comparing to a non-dynamic bin pool using Fisher???s exact test.
#The non-dynamic bin pool was sampled with replacement to
#match the distribution of average signal strength from the dynamic bins.
#Following that, significant TF PWMs were grouped in subfamilies using
#the structural information from TFClass71 because they share similar
#if not identical binding motifs. The top significantly over-represented TFs
#and their associated subfamilies were reported.

#Motif enrichment and GO analysis. Motif enrichment for each cell type.
# Motif enrichment for each cell type and histone modifications
# were carried out using ChromVAR61. Briefly, mapped reads were converted to
# cell-to-bin matrices with a bin size of 1,000 bp for four histone profiles.
# Reads for each bin were summarized from all cells of the same
# groups from transcriptome-based clustering. GC bias and background peaks were
# calculated and motif enrichment score for each cell type was
# then computed using the computeDeviations function of ChromVAR61.
# Motif enrichment for each CRE module: Motif enrichment for each CRE module
# was analyzed using Homer (v.4.11)62. We scanned a region of ??200 bp
# around the center of the element for both de novo and known motif enrichment analysis
# (from the JASPAR database63). We used the total peak list as the background
# for motif enrichment analysis of cCREs in each group.
#
#61.Schep, A. N., Wu, B., Buenrostro, J. D. & Greenleaf, W. J. chromVAR: inferring
#transcription-factor-associated accessibility from single-cell epigenomic data. Nat. Methods 14, 975???978 (2017).

#62.Heinz, S. et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory
#elements required for macrophage and B cell identities. Mol. Cell 38, 576???589 (2010).
