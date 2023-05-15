library("rtracklayer")
library("dbplyr")
library("dplyr")
library("motifStack")
library("universalmotif")
library(ggplot2)


#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"

gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds"))



all_motif_matches_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_Teemu_rest.Rds")
length(names(all_motif_matches_Teemu)) #15
length(unique(names(all_motif_matches_Teemu))) #15

motifs=names(all_motif_matches_Teemu) #15
match_numbers=unlist(lapply(all_motif_matches_Teemu, length))

saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu_rest.Rds")

#relaxed affinity threshold
motif_matches_4_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_4_Teemu.Rds") #1025
which(table(names(motif_matches_4_Teemu))!=1)
#For some reason, "ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES" is twice

motif_matches_4_Teemu=motif_matches_4_Teemu[-grep("ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES", names(motif_matches_4_Teemu))[2]]
#all_motif_matches_Teemu[grep("ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES", names(all_motif_matches_Teemu))] it is not twice in the original data

#

#There are some motifs which occur multiple times
motif_matches_3_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_3_Teemu.Rds") #762
length(names(motif_matches_3_Teemu)) #762
length(unique(names(motif_matches_3_Teemu))) #761
which(table(names(motif_matches_3_Teemu))!=1)
motif_matches_3_Teemu=motif_matches_3_Teemu[-grep("ELK1_PAX5_CAP-SELEX_TCAGTT40NTTG_AAA_RSCGGAACYACGCWYSANTG_2_2_YES", names(motif_matches_3_Teemu))[2]]


motifs=names(all_motif_matches_Teemu) #3892
motifs_4=names(motif_matches_4_Teemu) #1024
motifs_3=names(motif_matches_3_Teemu) #763
length(unique(motifs_3)) #761
which(table(motifs_3)>1)

match_numbers=unlist(lapply(all_motif_matches_Teemu, length))
match_numbers_4=unlist(lapply(motif_matches_4_Teemu, length))
match_numbers_3=unlist(lapply(motif_matches_3_Teemu, length))

#saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds")
#saveRDS(match_numbers_4,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_4_Teemu.Rds")
#saveRDS(match_numbers_3,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_3_Teemu.Rds")
match_numbers <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds")
match_numbers_4 <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu_4.Rds")
match_numbers_3 <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds")


# Basic histogram

df=data.frame(match_nro=match_numbers)
max(df$match_nro) #Vierstra: 4063189 Teemu: 758738
min(df$match_nro) #Vierstra: 15846 Teemu: 216

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu_motif-hit-nros.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=match_nro)) + geom_histogram(binwidth=5000)+ylim(c(0,300))
dev.off()

df=data.frame(match_nro=match_numbers_4)
max(df$match_nro) #Vierstra: 4063189 Teemu: 758738 Teemu_4 358533
min(df$match_nro) #Vierstra: 15846 Teemu: 216 Teemu_4 337

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu_4_motif-hit-nros.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=match_nro)) + geom_histogram(binwidth=5000)
dev.off()

df=data.frame(match_nro=match_numbers_3)
max(df$match_nro) #Vierstra: 4063189 Teemu: 758738 Teemu_4 358533 Teemu_3 304471
min(df$match_nro) #Vierstra: 15846 Teemu: 216 Teemu_4 337 Teemu_3 480

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu_3_motif-hit-nros.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=match_nro)) + geom_histogram(binwidth=20000)
dev.off()


#For which motifs the length of matches is less than 300 000


rerun=which(match_numbers < 300000) #1024  CTCF: 51136 
rerun4=which(match_numbers_4 < 300000) #761 CTCF: 77283
rerun3=which(match_numbers_3 < 300000) #578 CTCF 117589

hist(match_numbers[rerun])
hist(match_numbers[rerun4])
hist(match_numbers_3[rerun3])



write.table(names(rerun), "../../Experiments/rerunMOODS_lower_threshold.txt")
write.table(names(rerun4), "../../Experiments/rerunMOODS_lower_threshold_4.txt")
write.table(names(rerun3), "../../Experiments/rerunMOODS_lower_threshold_3.txt")

#hist(unlist(lapply(all_motif_matches_Vierstra, length)))





#How many of the rerun are representatives
#There are some duplicates in the motif matches
length(unique(names(rerun3))) #578
length(names(rerun3)) #580

which(table(names(rerun3))!=1)
#ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES 
#FLI1_CEBPB_CAP-SELEX_TTCGGG40NGAG_AAC_RNCGGANNTTGCGCAAN_1_3_NO 
#which(representatives$ID=="ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES") 
#which(representatives$ID=="FLI1_CEBPB_CAP-SELEX_TTCGGG40NGAG_AAC_RNCGGANNTTGCGCAAN_1_3_NO") 

rerun3[which(names(rerun3)=="ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES")]
rerun3[which(names(rerun3)=="FLI1_CEBPB_CAP-SELEX_TTCGGG40NGAG_AAC_RNCGGANNTTGCGCAAN_1_3_NO")]


length(which(names(rerun) %in% rep_motifs))
length(which(names(rerun4) %in% rep_motifs))
length(which(names(rerun3) %in% rep_motifs)) #182
length(which( rep_motifs %in% names(rerun3))) #181

head(representatives$ID)
representatives$match_numbers=match_numbers[match( test,names(match_numbers),)]
  
#head(names(match_numbers[match( test,names(match_numbers),)]))

representatives$enoughHits="YES"
representatives$enoughHits[representatives$match_numbers< 300000]="NO"

monomers=c("monomeric", "monomer or dimer", "monomer", "monomeric or dimeric")

multimers=c("dimeric", "putative multimer", "putatively multimeric", "dimer of dimers", "trimeric")

#both monomeric

representatives$motif_type="monomer"
cap_selex=representatives$experiment=="CAP-SELEX"
mm=representatives$type %in% multimers
representatives[mm | cap_selex, "motif_type"]="multimer"

library("ggplot2")

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-nros-monomer-multimer.pdf" ))
ggplot(representatives, aes(x=match_numbers, color=motif_type)) + geom_histogram(fill="white",binwidth=10000,aes(y=..density..))+ 
  xlim(0, 1e6)+
  theme_classic()
dev.off()

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-ICs.pdf" ))
ggplot(representatives, aes(x=IC, color=enoughHits)) + 
  geom_histogram(fill="white",binwidth=0.5, alpha=0.5, position="identity",aes(y=..density..)) +
scale_color_grey()+scale_fill_grey() +
  theme_classic()
dev.off()

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-lengths.pdf" ))
ggplot(representatives, aes(x=length, color=enoughHits)) + 
  stat_bin(aes(y=..density..), position='identity', fill="white", binwidth=1, alpha=0.5) +
  #geom_histogram(fill="white",binwidth=1, alpha=0.5, position="identity") +
  scale_color_grey()+scale_fill_grey() +
  theme_classic()
dev.off()

representatives$norm_IC=representatives$IC/representatives$length

#IC divided by the motif length
pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-ICs-normalized.pdf" ))
ggplot(representatives, aes(x=norm_IC, color=enoughHits)) + 
  stat_bin(fill="white", bins=50,alpha=0.5, position="identity",aes(y=..density..)) +
  scale_color_grey()+scale_fill_grey() +
  theme_classic()
dev.off()


ggplot(representatives, aes( x=IC,y=match_numbers))+geom_point()
ggplot(representatives, aes( x=norm_IC,y=match_numbers))+geom_point()

#scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")

#Extract the metadata for these motifs, list the number of motif matches with threshold 3 and motif lengths

rep_ind=which( rep_motifs %in%  names(rerun))

representatives_problem=representatives[which(representatives$new_representative=="YES"),]
representatives_problem=representatives_problem[rep_ind,]

#Are the motif lengths or information contents somehow different for motifs for which we obtain at least 300,000 matches 
#vs motifs for which we obtain less

#Add hit nro and motif length
representatives_problem$hit_nro=match_numbers[rep_motifs[rep_ind] ]
#representatives_problem$motif_length=sapply(motif_matches_3_Teemu[rep_motifs[rep_ind]], function(x){ unique(nchar(x$sequence)) })

write.table(representatives_problem, file="representatives_problem.tsv", sep="\t", row.names=FALSE, quote=FALSE)
representative_motif_matches_Teemu=all_motif_matches_Teemu[rep_motifs]
saveRDS(representative_motif_matches_Teemu, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")

representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")




rep_ind=which( rep_motifs %in%  names(rerun3))

#Add hit nro and motif length
representatives_problem$hit_nro=match_numbers_3[rep_motifs[rep_ind] ]
representatives_problem$motif_length=sapply(motif_matches_3_Teemu[rep_motifs[rep_ind]], function(x){ unique(nchar(x$sequence)) })

write.table(representatives_problem, file="representatives_problem.tsv", sep="\t", row.names=FALSE, quote=FALSE)
representative_motif_matches_Teemu=all_motif_matches_Teemu[rep_motifs]
saveRDS(representative_motif_matches_Teemu, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")

representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")

motifs=names(representative_motif_matches)


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


saveRDS(p_matrix,  file = paste0( "../RData/Teemu_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/Teemu_e_matrix_human_cell_type_restricted.Rds")) #





