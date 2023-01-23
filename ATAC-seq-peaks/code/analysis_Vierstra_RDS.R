library("rtracklayer")
library("dbplyr")
library("dplyr")

#Database stuff
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")


#all_motif_matches_Vierstra <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_Vierstra.Rds")

#motifs=names(all_motif_matches_Vierstra) #3892

#Which motifs are representative

#Which are representative motifs

#representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
#test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
#rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new

#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

#length(unique(rep_motifs))
#1062

#representative_motif_matches_Vierstra=all_motif_matches_Vierstra[rep_motifs]
#saveRDS(representative_motif_matches_Vierstra, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Vierstra.Rds")

representative_motif_matches_Vierstra <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Vierstra.Rds")

motifs=names(representative_motif_matches_Vierstra)


cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

p_matrix=matrix(nrow=length(motifs), ncol=length(cCREs_list),
                dimnames = list( motifs,
                                names(cCREs_list) )) #1062 x 111

e_matrix=matrix(nrow=length(motifs), ncol=length(cCREs_list),
                dimnames = list(motifs,
                                names(cCREs_list) )) #1062 x 111



all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs)


N=length(all_cCREs) #all possible cCREs #435142

#m: a list of all cCREs that have a motif match#

#select only seqnames, start, end from the table
#head(data.frame(rep_motifs, sort(as.vector(rep_id$PWMs)) ))
#table(rep_motifs==sort(as.vector(rep_id$PWMs)) ) all TRUE

#library(ggplot2)
# Basic histogram

#df=data.frame(cCREs_nro= unlist(lapply(cCREs_list, length)))

#max(df$cCREs_nro) #48391
#min(df$cCREs_nro) #846

#ggplot(df, aes(x=_nro)) + geom_histogram()
# Change the width of bins
#pdf(file = "../Figures/cCREs-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
#ggplot(df, aes(x=cCREs_nro)) + geom_histogram(binwidth=1000)
#dev.off()

#Width of cCREs, these are all the same 400 bp?

#widths=lapply(cCREs_list, width)
#lapply(widths, table)


for(hits in motifs){
  print(hits)
  #hits=motifs[1]

 
  #m: the size of all cCREs that have a motif match#
  fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  
  m=length(unique(fol@to)) #This is now number of all cCREs that have at least one motif hit #9857
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #findOverlaps(query, subject)
    fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
    #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
    # number of cCREs that have motif match
    
    k=length(unique(fol@to)) #This is now number of cCREs that have at least one motif hit #178
    #k=length(unique(fol@from)) This is now number of hits that overlap with cCREs (there can be multiple motifs hits in a single CRE)
    
    #m	the number of white balls in the urn.
    # number of motif matches, 300 000
    
    # N-m
    # the number of black balls in the urn.
    # the number of all possible cCREs - the number of motif matches
    
    #Nm=N-m #list of all CREs that do not have a motif This is computed in FishersExact
    
    # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
    # The number of cell type specific cCREs
    
    n=length(cCREs_list[[cell_type]]) #This is correct
    
    #Hypergeometric test/Fisher's exact test
    
    #m: a list of all cCREs that have a motif match
    
    #n: A set of cCREs restrictive to a given cell type. G_0, n=|G_0|
    
    #N: G: all cCREs across all cell types G, N=|G|
    
    #k: number of cell-line-specific cCREs that have motif match
    
    
    fe=FishersExact(k,n,m,N)
    p_matrix[hits, cell_type]=fe$p
    e_matrix[hits, cell_type]=fe$e
    
    
    
  }
  
}


saveRDS(p_matrix,  file = paste0( "../RData/Vierstra_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/Vierstra_e_matrix_human_cell_type_restricted.Rds")) #





