invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
library(rtracklayer)
library(furrr)

options(future.availableCores.methods = "Slurm")

availableCores()

# With plan(multicore) the number of workers is based on batch job reservation details.
plan("multicore")


setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
source("../code/FishersExact.R")

cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

GR_list_representative <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches.Rds")

#save smaller representative list for testing purposes
GR_list_small=GR_list_representative[1:4]
saveRDS(CREmodules,  file ="../RData/CREmodules.Rds") 

m=length(GR_list_representative[[hits]])

funtorun <- function(hits, cCREs_list){ 
  
  m=length(GR_list_representative[[hits]])
  
  p_vector[hits, cell_type]=fe$p
  e_matrix[hits, cell_type]=fe$e
  
  ?vector
  p_vector=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                  dimnames = list(names(GR_list_representative),
                                  names(cCREs_list) )) #1062 x 111
  
  e_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                  dimnames = list(names(GR_list_representative),
                                  names(cCREs_list) )) #1062 x 111
  
  
  
  
  for(cell_type in names(cCREs_list)){
  #print(cell_type)
  fol=findOverlaps(hits, cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
  k=length(unique(fol@to))
  Nm=N-m
  n=length(cCREs_list[[cell_type]])
  fe=FishersExact(k,n,m,N)
  
}






#compute the rows of these matrices in parallel
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
