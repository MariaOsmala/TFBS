library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"

gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds"))

#Motif metadata with ICs and length

representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length.tsv", sep="\t", header=TRUE)

tmp=c("NR1I2_HT-SELEX_TACTGG40NGGA_KR_AGTTCRNNNRGTTCA_1_4b2_NA",  "NR1I3_HT-SELEX_TTTAGC40NGAC_KL_AGTTCANNNNAGTTCA_1_3b1_NA",
 "ONECUT1_HT-SELEX_TAGCTC20NTCT_Y_CRATCRATAWN_1_3_NA",       "ONECUT2_HT-SELEX_TGGGCG30NCGT_AH_NNCGATCRATAWNN_1_2_NA"  ,
 "PAX3_HT-SELEX_TTAGGG20NGGA_AL_ASCGTGACTAATT_1_3_NA" ,      "PAX3_HT-SELEX_TTAGGG20NGGA_AL_NNTAATYGATTANN_1_3_NA"     ,
 "PAX4_HT-SELEX_TCGGGG20NGA_KG_NKAATTAGNTAATTAN_2_3_NA",     "PAX4_HT-SELEX_TCGGGG20NGA_KG_SGTCACGCATSRNTG_2_3_NA"     ,
 "PGR_HT-SELEX_TCGCAA40NTGT_AR_RGNACRNNNTGTNCY_1_4_NA",      "SREBF2_HT-SELEX_TGAGAT20NGA_Y_NRTCACGCCAYN_1_3_NA")

representatives$ID[which(representatives$ID %in% tmp )]=gsub("_NA", "",representatives$ID[which(representatives$ID %in% tmp )])

#representatives$ID[7]
#"ATF2_FLI1_TACAGA40NAGAT_YYIII_NTGACGTMACCGGAWRYN_m2_c3b0_short_pfm_composite_new"

#representatives$filename[7]
#"PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/ATF2_FLI1_TACAGA40NAGAT_YYIII_NTGACGTMACCGGAWRYN_m2_c3b0_short.pfm"

#filenames[1]
#"ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short_composite.pfm"

#rename Yimengs motifs, add 

#filenames=system("ls ../../PWMs/fromYimeng/pwms_space/pfm_composite_new/*.pfm")
filenames=dir("../../PWMs/fromYimeng/pwms_space/pfm_composite_new/")
match_with_ID=gsub("short_composite", "short_pfm_composite_new", filenames)
match_with_ID=gsub("\\.pfm", "", match_with_ID)

match_ind=match( match_with_ID,representatives$ID)


representatives$ID[match_ind]=gsub(".pfm", "", filenames)

#Maybe check the filenames at some point
#representatives$filename[match_ind]="../../PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/" filenames

filenames2=dir("../../PWMs/fromYimeng/pwms_space/pfm_spacing_new/")

match_with_ID=gsub("short_spacing", "short_pfm_spacing_new", filenames2)
match_with_ID=gsub("\\.pfm", "", match_with_ID)

match_ind=match( match_with_ID,representatives$ID)
representatives$ID[match_ind]=gsub(".pfm", "", filenames2)
#write.table(representatives, 
#file="/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames.tsv", 
#      quote=FALSE, sep="\t", row.names=FALSE)
#test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])


rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

#grep("CTCF",rep_motifs) CTCF is in representative motif

#remove _pfm_composite_new and _pfm_spacing_new
#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)
#length(unique(rep_motifs))
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




#load( file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_5_4_3_2_Teemu.RData") #1025
#all_motif_matches_Teemu,match_numbers,

#match_numbers=unlist(lapply(all_motif_matches_Teemu, length))

#rep_motifs[which( !(rep_motifs %in% names(all_motif_matches_Teemu)))]=gsub("_NA", "", rep_motifs[which( !(rep_motifs %in% names(all_motif_matches_Teemu)))] )

#representative_motif_matches_Teemu=all_motif_matches_Teemu[rep_motifs]

#saveRDS(representative_motif_matches_Teemu, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")



representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Teemu.Rds")
#match_numbers=unlist(lapply(representative_motif_matches, length))

pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/Teemu-match-number-distribution-representatives.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
df=data.frame(match_nro=match_numbers)
ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_cut(breaks=c(0.03), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, top 300000 with threshold 2, binwidth 20000')+ theme_minimal()
dev.off()


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

saveRDS(p_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/Teemu_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/Teemu_e_matrix_human_cell_type_restricted.Rds")) #





