library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
#.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))
library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"
#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

#grep("CTCF",rep_motifs) CTCF is in representative motif

length(unique(rep_motifs))
#1031

#top_motif_matches_human_final=readRDS(file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds") # 3124? should be 3294, more?


#head(names(top_motif_matches_human_final)[which(!(names(top_motif_matches_human_final) %in% representatives$ID))])

#missing=representatives$ID[which(!(representatives$ID %in% names(top_motif_matches_human_final) ))] #170

#all_motif_matches_Teemu,match_numbers,

#match_numbers=unlist(lapply(top_motif_matches_human_final, length))
#saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_top_motif_matches_human_final.Rds") 

#match_numbers=readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_top_motif_matches_human_final.Rds")


#df=data.frame(ID=names(match_numbers),match_numbers=as.vector(match_numbers))

# Use the merge function to add a column from df2 to df1
#merged_representatives <- merge(representatives, df[,c('ID','match_numbers')], by = 'ID', all.x = TRUE)

#Minimum integer threshold used to call the motifs

# floor(min(top_motif_matches_human_final[[1]]$score))
# f1 <- function(x) floor(min(x$score))
# thresholds=lapply(top_motif_matches_human_final, f1) #, simplity="array"
# 
# df=data.frame(ID=names(top_motif_matches_human_final),threshold=unlist(thresholds))
# 
# merged_representatives2 <- merge(merged_representatives, df[,c('ID','threshold')], by = 'ID', all.x = TRUE)
# 
# write.table(merged_representatives2, "../../PWMs_final/metadata_representatives_match_numbers_thresholds.tsv",sep="\t", 
#             col.names=TRUE, row.names=FALSE) #3294
# 



#Motif number distribution
# pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/motif-match-number-distribution-final.pdf" ))
# #width, height, onefile, family, title, fonts, version,
# #paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
# #useDingbats, useKerning, fillOddEven, compress)
# df=data.frame(match_nro=match_numbers)
# ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
#   scale_y_cut(breaks=c(0.05), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
#   #scale_y_break(c(0.0225, 0.93 ))
#   #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
#   labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, binwidth 20000')+ theme_minimal()
# dev.off()





#representative_motif_matches=top_motif_matches_human_final[rep_motifs]
#saveRDS(representative_motif_matches, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")
representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")
# match_numbers_representatives=unlist(lapply(representative_motif_matches, length))
# 
# pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/motif-match-number-distribution-final-representatives.pdf" ))
# #width, height, onefile, family, title, fonts, version,
# #paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
# #useDingbats, useKerning, fillOddEven, compress)
# df=data.frame(match_nro=match_numbers_representatives)
# ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
#   scale_y_cut(breaks=c(0.03), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
#   #scale_y_break(c(0.0225, 0.93 ))
#   #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
#   labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, top 300000 with threshold 2, binwidth 20000')+ theme_minimal()
# dev.off()


#match_numbers[grep("CTCF", names(match_numbers))] #177144


motifs=names(representative_motif_matches)


cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
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

for(hits in motifs){
  print(hits)
  #hits=motifs[1]

 
  #m: the size of all cCREs that have a motif match#
  fol=findOverlaps(representative_motif_matches[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  
  m=length(unique(fol@to)) #This is now number of all cCREs that have at least one motif hit #9857
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #findOverlaps(query, subject)
    fol=findOverlaps(representative_motif_matches[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
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

saveRDS(p_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted.Rds")) #





