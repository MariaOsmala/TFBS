library("rtracklayer")
library("dbplyr")
library("dplyr")


setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
source("../code/FishersExact.R")
gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")


# * Cell-type specificity of cCREs across the cell types.
# 1.2 million cCREs clustered into 150 cis-regulatory modules (CRMs).
# Putative master regulators of fetal and adult cell types.

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
# library(readr)
# gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")
# 
# cCREs=read_delim("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_cCRE_list_GRCh38.tsv.gz", delim="\t", col_names=TRUE)
# nrow(cCREs) #1154611, this is the same table as cCRE_hg38
# 
# names(cCREs)=c("chr", "start", "end", "class", "fetal", "adult", "CREmoduleID")
# 
# library("tidyverse")
# library("GenomicRanges")
# modules=cCREs %>%
#   group("CREmoduleID")
# 
# by_CREdmodules <- cCREs %>% group_by(CREmoduleID)
# 
# by_CREdmodules %>% tally(sort=TRUE)
# 
# group_indices=by_CREdmodules %>% group_indices()
# 
# CREmodules <- GRangesList()
# 
# for(gk in 1:150){
#   print(gk)
#   
#   data=as.data.frame(cCREs %>% filter( CREmoduleID==gk) %>% select(c("chr", "start", "end")))
#   tmp_GRanges=GRanges(seqnames = Rle( data$chr, rep(1, nrow(data)) ),
#                       ranges = IRanges(start=data$start, end = data$end) ,
#                       strand = Rle(strand( rep("*", nrow(data))) , rep(1, nrow(data))  ))
#   
#   
#   isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
#   
#   seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
#   names(seqlengths)=paste0("chr", names(seqlengths))
#   
#   ############## try dropSeqlevels or ?keepStandardChromosomes from GenomicRanges library.
#   
#   
#   tmp_GRanges=dropSeqlevels(tmp_GRanges, names(seqlengths(tmp_GRanges)[ which( (as.vector(names(seqlengths(tmp_GRanges))) %in% names(seqlengths)==FALSE) ) ]),
#                             pruning.mode = "coarse")
#   
# 
#   
#   seqlengths(tmp_GRanges)=seqlengths[names( seqlengths(tmp_GRanges))]
#   
#   genome(tmp_GRanges) <- "GRCh38"
#   
#   CREmodules[[gk]]=tmp_GRanges
#   
#   
# }

#saveRDS(CREmodules,  file ="../RData/CREmodules.Rds") 

CREmodules<-readRDS(file = "../RData/CREmodules.Rds")

#Modules connected to cell_types, groups of cell types

# Meta_modules=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5D_CRE_Module_Motif_enrichment.csv",sep="\t")
# Meta_modules=as.data.frame(t(Meta_modules[c(1,2), -1]))
# rownames(Meta_modules)=NULL
# names(Meta_modules)=c("ModuleGroup", "ModuleID")
# Meta_modules[,"ModuleID"]=as.numeric(Meta_modules[,"ModuleID"])
# 
# saveRDS(Meta_modules,  file ="../RData/Meta_modules.Rds") 

con <- DBI::dbConnect(RSQLite::SQLite(), "../SQLite/motif-hits.sqlite")

representative <- tbl(con, "representative")

#Unique motifs in the representative table
motif_names <- representative %>%  distinct(PWM) %>% collect()

p_matrix=matrix(nrow=length(motif_names$PWM), ncol=length(CREmodules),
                dimnames = list( motif_names$PWM,
                                 1:length(CREmodules) )) #1062 x 111

e_matrix=matrix(nrow=length(motif_names$PWM), ncol=length(CREmodules),
                dimnames = list(motif_names$PWM,
                                1:length(CREmodules) )) #1062 x 111

N=length(unique(unlist(CREmodules))) #all possible cCREs
# 1_154_464

for(hits in as.vector(motif_names$PWM)){
  print(hits)
  #hits=as.vector(motif_names$PWM)[1]
  
  
  res=DBI::dbSendQuery(con, paste0("SELECT seqnames, start, end FROM `representative` WHERE (`PWM` IN ('",hits,"'))"))
  data=DBI::dbFetch(res)
  
  DBI::dbClearResult(res)
  
  tmp_GRanges=GRanges(seqnames = Rle( data$seqnames, rep(1, nrow(data)) ),
                      ranges = IRanges(start=data$start, end = data$end) ,
                      strand = Rle(strand( rep("*", nrow(data))) , rep(1, nrow(data))  ))
  
  
  isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  
  seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  
  seqlengths(tmp_GRanges)=seqlengths[names( seqlengths(tmp_GRanges))]
  genome(tmp_GRanges) <- "GRCh38"
  
  #m	the number of white balls in the urn.
  # number of motif matches, 300 000
  m=length(tmp_GRanges)
  
  for(cell_type in 1:length(CREmodules)){
    #print(cell_type)
    
    
    fol=findOverlaps(tmp_GRanges, CREmodules[[cell_type]], type="within", ignore.strand=TRUE)
    #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
    # number of cCREs that have motif match
    
    k=length(unique(fol@to))
    
  
    
    # N-m
    # the number of black balls in the urn.
    # the number of all possible cCREs - the number of motif matches
    
    Nm=N-m
    
    # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
    # The number of cell type specific cCREs
    
    n=length(CREmodules[[cell_type]])
    
    fe=FishersExact(k,n,m,N)
    p_matrix[hits, cell_type]=fe$p
    e_matrix[hits, cell_type]=fe$e
    }
  
}


DBI::dbDisconnect(con)
saveRDS(p_matrix,  file = paste0( "../RData/p_matrix_CRE_modules.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/e_matrix_CRE_modules.Rds")) #






# * Phylogenetic analysis to group fetal and adult cell types into different cell lineages, enrichment for lineage-specific TFBSs?
#   * Differentially accessible (DA) elements between adults and fetal cell lines/cell lineage groups.
#    Adult and fetal-specific cCREs for the cell-lineage groups. Specific cCREs underlying fetal or adult-specific regulatory programs.

#WHERE IS THIS DATA?

