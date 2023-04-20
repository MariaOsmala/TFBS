library(rtracklayer)

#Database stuff
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
source("../code/FishersExact.R")


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
# 7. CRE module: The ID of CRE module that the cCRE belongs to.

# * 435,142 adult cCREs with restricted accessibility in one or a few cell types

#From Mendeley database
#files=dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/")

#ct_names=do.call(rbind,sapply(files, strsplit, ".bed.gz"))
#dimnames(ct_names)=NULL

#cCREs_list<-GRangesList()

#for(ct in ct_names){
#  print(ct)
#  cCREs_list[[ct]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/",ct,".bed.gz"),format="bed")


#}

#saveRDS(cCREs_list,  file = paste0( "../RData/cCREs_list.Rds")) #

#sum(sapply(cCREs_list, length)) #1,091,805

#unique
#length(unique(unlist(cCREs_list))) #435142

cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

#Which are representative motifs

#representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
#test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
#rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and
#_pdf_spacing_new

#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

#length(unique(rep_motifs))
#1062
#setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject")
#GR_list<-readRDS(file = paste0( "../RData/all_motif_matches_sbatch.Rds")) #
#length(GR_list)
#3982

#GR_list_representative=GR_list[ rep_motifs]
#GR_list_representative=GR_list[ which ( ( names(GR_list) %in% rep_motifs)==TRUE)]

#saveRDS(GR_list_representative,  file = paste0( "../RData/representative_motif_matches.Rds")) #


GR_list_representative <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches.Rds")



p_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                dimnames = list(names(GR_list_representative),
                                names(cCREs_list) )) #1062 x 111

e_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                dimnames = list(names(GR_list_representative),
                                names(cCREs_list) )) #1062 x 111

N=length(unique(unlist(cCREs_list))) #all possible cCREs

for(hits in names(GR_list_representative)){
    print(hits)
    #hits=names(GR_list_representative)[1]
    m=length(GR_list_representative[[hits]])
    for(cell_type in names(cCREs_list)){
        print(cell_type)
        #cell_type=names(cCREs_list)[1]
        #findOverlaps(query=hits, subjects=cCREs_list[cell_type])
      
        #If type is within, the query interval must be wholly contained within 
        #the subject interval. Note that all matches must additionally satisfy 
        #the minoverlap (default 0) constraint described above.
      
        #The maxgap parameter has special meaning with the special overlap types. 
        #For within, it is the maximum amount by which the subject may be wider than the query. 
        #If maxgap is set to -1 (the default), it's replaced internally by 0.
      
      
        fol=findOverlaps(GR_list_representative[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
        #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
        # number of cCREs that have motif match
    
        k=length(unique(fol@to))
    
        #m	the number of white balls in the urn.
        # number of motif matches, 300 000
        
    
    
    
        # N-m
        # the number of black balls in the urn.
        # the number of all possible cCREs - the number of motif matches
    
        Nm=N-m
    
        # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
        # The number of cell type specific cCREs
    
        n=length(cCREs_list[[cell_type]])
        
        fe=FishersExact(k,n,m,N)
        p_matrix[hits, cell_type]=fe$p
        e_matrix[hits, cell_type]=fe$e
        
    }
  
}
  
saveRDS(p_matrix,  file = paste0( "../RData/p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/e_matrix_human_cell_type_restricted.Rds")) #
  
  
  


# * Phylogenetic analysis to group fetal and adult cell types into different cell lineages, enrichment for lineage-specific TFBSs?
#   * Differentially accessible (DA) elements between adults and fetal cell lines/cell lineage groups.
#    Adult and fetal-specific cCREs for the cell-lineage groups. Specific cCREs underlying fetal or adult-specific regulatory programs.




# * Cell-type specificity of cCREs across the cell types.
# 1.2 million cCREs clustered into 150 cis-regulatory modules (CRMs).
# Putative master regulators of fetal and adult cell types.



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
# library("Matrix")
# matrix=readMM("../CATLAS/matrix.tsv.gz") #this is 1,154,611 x 222 matrix, contains TRUE or FALSE
# 
# tmp=do.call(rbind, apply(matrix,1,table)) #For each cCREs, in how many it is present
# restricted_cCREs=which(tmp[,"TRUE"]==1) #cCREs that are present in a single cell type  347991 (paper: 345142 cCREs, maybe blacklist removed)
# 
# celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t")

#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#library(rtracklayer)



#cCREs=import("../CATLAS/cCREs.bed.gz",format="bed")





#read TF peaks
#Yimengs motifs
#RFX3 GLI3
#FLI1 GLI3
#TCF7

#TCF7=read.table("../../Results/MOODS/TCF7.csv",sep=",")

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
