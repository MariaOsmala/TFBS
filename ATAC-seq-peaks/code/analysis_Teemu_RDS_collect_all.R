library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library(ggbreak)
#install.packages("ggbreak")



setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"

gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds"))

#Motif metadata with ICs and length

#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length.tsv", sep="\t", header=TRUE)
representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames.tsv", sep="\t", header=TRUE)


test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new
rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)
length(unique(rep_motifs))
#1062

#Interesting motifs for Yimen's paper, choose the representative ones, missing?

#interesting=c("PTF1A_NKX6.1",
#"OLIG2_OTP",
#"OLIG2_EVX2",
#"NEUROG2_EVX2",
#"ATOH1_EVX2",
# "BHLHE22_EVX2",
# "NHLH1_EVX2",
# "OLIG2_HOXB5",
# "OLIG2_NKX6.1",
# "ATOH1_EVX2",
# "NHLH1_HOXD8",
# "PROX1_HOXB2")

# representatives$symbol[grep("PTF1A", representatives$symbol)]
# representatives$symbol[grep("NKX6", representatives$symbol)]
# representatives$symbol[grep("OLIG2", representatives$symbol)]
# representatives$symbol[grep("OTP", representatives$symbol)]
# representatives$symbol[grep("NEUROG2", representatives$symbol)]
# representatives$symbol[grep("EVX2", representatives$symbol)]
# representatives$symbol[grep("ATOH1", representatives$symbol)]
# representatives$symbol[grep("BHLHE22", representatives$symbol)]
# representatives$symbol[grep("NHLH1", representatives$symbol)]
# representatives$symbol[grep("OLIG2", representatives$symbol)]
# representatives$symbol[grep("HOXB5", representatives$symbol)]
# representatives$symbol[grep("ATOH1", representatives$symbol)]
# representatives$symbol[grep("NHLH1", representatives$symbol)]
# representatives$symbol[grep("HOXD8", representatives$symbol)]
# representatives$symbol[grep("PROX1_HOXB2", representatives$symbol)]
# 
# interesting_df=representatives[which(representatives$symbol %in% interesting),]


all_motif_matches_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_Teemu.Rds") #50 GB

length(names(all_motif_matches_Teemu)) #3982
length(unique(names(all_motif_matches_Teemu))) #3982
#match_numbers=unlist(lapply(all_motif_matches_Teemu, length)) CHECK!
#saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds") CHECK!
match_numbers <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds")

df=data.frame(match_nro=match_numbers)
max(df$match_nro) # Teemu: 758738
min(df$match_nro) # Teemu: 216

#THese are only the top, that's why sharp break around 300000
pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu_motif-hit-nros.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=match_nro)) + geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_continuous(labels = scales::percent)+ 
  scale_y_cut(breaks=c(0.1), which=c(1,2), scales=c(0.5,0.5))+ylim(0,0.8)+#xlim(c(280000,325000))
  scale_x_continuous(labels = scales::comma)+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers with threshold 5, binwidth 20000')+ theme_minimal()
dev.off()


#Remove those with not enought matches

rerun=which(match_numbers < 300000) #1024  CTCF: 51136

all_motif_matches_Teemu=all_motif_matches_Teemu[-rerun]

match_numbers=match_numbers[-rerun]

#Add info to representatives

representatives$match_numbers=NA
ind=match( names(match_numbers),representatives$ID,) #3982-1024=2958

head(match_numbers[ind])
representatives$match_numbers[ind]=as.numeric(match_numbers)

representatives$MOODS_threshold=NA
representatives$MOODS_threshold[ind]=5

#relaxed affinity threshold
motif_matches_4_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_4_Teemu.Rds") #1025
which(table(names(motif_matches_4_Teemu))!=1)
#For some reason, "ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES" is twice

motif_matches_4_Teemu=motif_matches_4_Teemu[-grep("ETS1_HT-SELEX_TGTACC30NCAG_AI_NCCGGAWRYRYWTCCGGN_1_3_YES", names(motif_matches_4_Teemu))[2]]

#match_numbers_4=unlist(lapply(motif_matches_4_Teemu, length)) CHECK!
#saveRDS(match_numbers_4,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_4_Teemu.Rds") CHECK!
match_numbers_4 <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_4_Teemu.Rds")

all_motif_matches_Teemu=c(all_motif_matches_Teemu, motif_matches_4_Teemu)
length(all_motif_matches_Teemu)

match_numbers=c(match_numbers, match_numbers_4)

save(all_motif_matches_Teemu,match_numbers, file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_Teemu.RData") #1025 #CHECK
#load("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_Teemu.RData") #1025

#update representatives table
rerun=which(match_numbers_4 < 300000) #761
match_numbers_4=match_numbers_4[-rerun]

ind=match( names(match_numbers_4),representatives$ID,) #3982-1024=2958

head(match_numbers_4)
head(representatives$ID[ind])
representatives$match_numbers[ind]=as.numeric(match_numbers_4)


representatives$MOODS_threshold[ind]=4


#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
df=data.frame(match_nro=match_numbers)
max(df$match_nro) # Teemu: 758738
min(df$match_nro) # Teemu: 337

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-nros-added-4.pdf" ))
ggplot(df, aes(x=match_nro)) + geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_continuous(labels = scales::percent)+ 
  scale_y_cut(breaks=c(0.1), which=c(1,2), scales=c(0.5,0.5)) + ylim(0, 0.8) + #xlim(c(280000,325000)) +
  scale_x_continuous(labels = scales::comma)+theme_minimal()+
#scale_y_break(c(0.0225, 0.93 ))
#geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers with threshold 5, 4, binwidth 20000')
dev.off()



motif_matches_3_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_3_Teemu.Rds") #762
length(names(motif_matches_3_Teemu)) #762
length(unique(names(motif_matches_3_Teemu))) #761
which(table(names(motif_matches_3_Teemu))!=1)

motif_matches_3_Teemu=motif_matches_3_Teemu[-grep("ELK1_PAX5_CAP-SELEX_TCAGTT40NTTG_AAA_RSCGGAACYACGCWYSANTG_2_2_YES", names(motif_matches_3_Teemu))[2]]

#match_numbers_3=unlist(lapply(motif_matches_3_Teemu, length)) CHECK!
#saveRDS(match_numbers_3,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_3_Teemu.Rds") CHECK!
match_numbers_3 <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Teemu.Rds")

rerun=which(match_numbers < 300000) #761


tmp=all_motif_matches_Teemu[-rerun]
all_motif_matches_Teemu=c(tmp, motif_matches_3_Teemu)

tmp_match_numbers=match_numbers[-rerun]
match_numbers=c(tmp_match_numbers, match_numbers_3)

save(all_motif_matches_Teemu,match_numbers, file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_3_Teemu.RData") #1025 #CHECK

#load("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_3_Teemu.RData")

rm(tmp)
rm(tmp_match_numbers)

df=data.frame(match_nro=match_numbers)


max(df$match_nro) # Teemu: 758738
min(df$match_nro) # Teemu: 480

match_numbers[grep("CTCF", names(match_numbers))] #117589


#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-nros-added-3.pdf" ))
ggplot(df, aes(x=match_nro)) + geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_continuous(labels = scales::percent)+ 
  scale_y_cut(breaks=c(0.1), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 0.8)) + #xlim(c(280000,325000))
  scale_x_continuous(labels = scales::comma)+theme_minimal()+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers with threshold 5, 4, 3, binwidth 20000')
dev.off()

rerun=which(match_numbers_3 < 300000) #578
match_numbers_3=match_numbers_3[-rerun]

ind=match( names(match_numbers_3),representatives$ID,) #

head(match_numbers_3)
head(representatives$ID[ind])
representatives$match_numbers[ind]=as.numeric(match_numbers_3)


representatives$MOODS_threshold[ind]=3


motif_matches_2_Teemu <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_2_Teemu.Rds") #762
length(names(motif_matches_2_Teemu)) #579
length(unique(names(motif_matches_2_Teemu))) #578
which(table(names(motif_matches_2_Teemu))!=1)


motif_matches_2_Teemu=motif_matches_2_Teemu[-grep("RFX3_SREBF2_CAP-SELEX_TAAAGG40NGAC_AY_NNRGYAACNTCACGTGAY_1_3_YES", 
                                                  names(motif_matches_2_Teemu))[2]]

#match_numbers_2=unlist(lapply(motif_matches_2_Teemu, length)) CHECK!
#saveRDS(match_numbers_2,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_2_Teemu.Rds") CHECK!

match_numbers_2=readRDS(match_numbers_2,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_2_Teemu.Rds")


rerun=which(match_numbers < 300000) #576
tmp=all_motif_matches_Teemu

all_motif_matches_Teemu=all_motif_matches_Teemu[-rerun]



all_motif_matches_Teemu=c(all_motif_matches_Teemu, motif_matches_2_Teemu)
length(unique(names(all_motif_matches_Teemu)))
rm(tmp)

#update match_numbers
match_numbers=match_numbers[-rerun]
match_numbers=c(match_numbers, match_numbers_2)

match_numbers[grep("CTCF", names(match_numbers))] #177144

save(all_motif_matches_Teemu,match_numbers, file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_3_2_Teemu.RData") #1025

#load("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_3_2_Teemu.RData")


ind=match( names(match_numbers_2),representatives$ID) #

head(match_numbers_2)
head(representatives$ID[ind])

representatives$match_numbers[ind]=as.numeric(match_numbers_2)
representatives$MOODS_threshold[ind]=2

df=data.frame(match_nro=match_numbers)

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-motif-hit-nros-added-2.pdf" ), width=10, height=7)
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)

ggplot(df, aes(x=match_nro)) + geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_continuous(labels = scales::percent)+ 
  scale_y_cut(breaks=c(0.1), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 0.8))+#xlim(c(280000,325000))
  scale_x_continuous(labels = scales::comma)+theme_minimal()+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers with threshold 5, 4, 3, 2, binwidth 20000')
dev.off()

#What is the super hitter

head(representatives[order(representatives$match_numbers,decreasing = FALSE), c("ID", "match_numbers", "length", "IC")],20)

write.table(representatives, 
            file="/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames_match_numbers_MOODS_threshold.tsv", 
                  quote=FALSE, sep="\t", row.names=FALSE)


rerun=which(match_numbers < 300000) #454

representatives$match_numbers=match_numbers

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
ggplot(representatives, aes(x=match_numbers, color=motif_type)) + 
  geom_histogram(fill="white",binwidth=20000,aes(y=after_stat(count / sum(count))))+ 
  scale_y_continuous(labels = scales::percent)+
  xlim(0, 1e6)+
  scale_y_cut(breaks=c(0.05), which=c(1,2), scales=c(0.5,0.5)) + #ylim(c(0, 0.9))+
  theme_classic()+
labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers')
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

rerun=which(representatives$match_numbers< 300000)

rep_ind=which( rep_motifs %in%  names(rerun))

representatives_problem=representatives[which(representatives$new_representative=="YES"),]
representatives_problem=representatives_problem[rep_ind,]

#Now there are no representative in those which do not have enought hits

ggplot( (representatives %>% filter(new_representative=="YES") %>% select(match_numbers)), aes(x=match_numbers))+geom_histogram(binwidth = 20000)

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





