# library("rtracklayer")
# library("dbplyr")
# library("dplyr")
#library("motifStack")
#library("universalmotif")
 library("ggplot2")
 library("readr")
.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")

library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

args <- commandArgs(trailingOnly = TRUE)
genome=args[1]
#genome="hg38.analysisSset"
#genome="T2T-CHM13"

data_path="/scratch/project_2006203/TFBS/"

#Motif metadata with ICs and length
representatives=read.table("/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/metadata_final.tsv",sep="\t", header=TRUE) #3933
representatives=read.table(paste0("/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/metadata_final_match_numbers_thresholds_",
                                  "hg38.analysisSset",".tsv"),sep="\t", header=TRUE) 

representatives=read.table("/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/metadata_final_match_numbers_thresholds_T2T-CHM13.tsv",
                           ,sep="\t", header=TRUE) 

rep_motifs=representatives$ID[which(representatives$representative=="YES")]

length(unique(rep_motifs))
#1232

top_motif_matches_human_final=readRDS(file=paste0("/scratch/project_2006203/TFBS/Results/RData/",genome,".Rds")) # 
print(length(top_motif_matches_human_final)) #3933




match_numbers=unlist(lapply(top_motif_matches_human_final, length))
saveRDS(match_numbers,paste0("/scratch/project_2006203/TFBS/Results/RData/match_numbers_",genome,".Rds")) 
#match_numbers=readRDS(paste0("/scratch/project_2006203/TFBS/Results/RData/match_numbers_",genome,".Rds"))

#df=data.frame(ID=names(match_numbers),match_number_hg38.analysis.Set=as.vector(match_numbers))
df=data.frame(ID=names(match_numbers),match_number_T2T_CHM13=as.vector(match_numbers))

# Use the merge function to add a column from df2 to df1
#merged_representatives <- merge(representatives, df[,c('ID','match_number_hg38.analysis.Set')], by = 'ID', all.x = TRUE)
merged_representatives <- merge(representatives, df[,c('ID','match_number_T2T_CHM13')], by = 'ID', all.x = TRUE)

#Minimum integer threshold used to call the motifs

floor(min(top_motif_matches_human_final[[1]]$score))
min(top_motif_matches_human_final[[1]]$score)

f1 <- function(x) floor(min(x$score))

thresholds=lapply(top_motif_matches_human_final, f1) #, simplity="array"
thresholds_double=lapply(top_motif_matches_human_final, function(x) min(x$score)) #, simplity="array"
# 
#df=data.frame(ID=names(top_motif_matches_human_final),
#              threshold_hg38.analysis.Set=unlist(thresholds), threshold_double_hg38.analysis.Set=unlist(thresholds_double))

df=data.frame(ID=names(top_motif_matches_human_final),
              threshold_T2T_CHM13=unlist(thresholds), threshold_double_T2T_CHM13=unlist(thresholds_double))
# 
# merged_representatives2 <- merge(merged_representatives, df[,c('ID','threshold_hg38.analysis.Set', 'threshold_double_hg38.analysis.Set')], 
#                                  by = 'ID', all.x = TRUE)
 
merged_representatives2 <- merge(merged_representatives, df[,c('ID','threshold_T2T_CHM13', 'threshold_double_T2T_CHM13')], 
                                  by = 'ID', all.x = TRUE) 
 
# 
write_tsv(merged_representatives2,
          paste0("/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/metadata_final_match_numbers_thresholds_",genome,".tsv")) 
# 

representative_motif_matches=top_motif_matches_human_final[rep_motifs]
#saveRDS(representative_motif_matches, "/scratch/project_2006203/TFBS/Results/RData/representative_motif_matches_hg38.analysis.Set.Rds")
saveRDS(representative_motif_matches, "/scratch/project_2006203/TFBS/Results/RData/representative_motif_matches_T2T_CHM13.Rds")
#Motif number distribution
pdf(file = paste0(data_path, "Figures/motif-match-number-distribution-hg38-repeat-masked.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
df=data.frame(match_nro=merged_representatives2$match_number)
ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_cut(breaks=c(0.05), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, binwidth 20000')+ theme_minimal()
dev.off()



#THis figure failed
pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/motif-match-number-distribution-final-representatives_version2.2.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
df=data.frame(match_nro=match_numbers_representatives)
ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_cut(breaks=c(0.03), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, top 300000 with threshold 2, binwidth 20000')+ theme_minimal()
dev.off()


