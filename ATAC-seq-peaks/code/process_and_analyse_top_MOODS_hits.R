library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library("readr")
.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")

library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")



data_path="/scratch/project_2006203/TFBS/"

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3933
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

#grep("CTCF",rep_motifs) CTCF is in representative motif

length(unique(rep_motifs))
#1232

top_motif_matches_human_final=readRDS(file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds") # 

print(length(top_motif_matches_human_final))

print(head(names(top_motif_matches_human_final)[which(!(names(top_motif_matches_human_final) %in% representatives$ID))]))


#match_numbers=unlist(lapply(top_motif_matches_human_final, length))
#saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_top_motif_matches_human_final_version2.2.Rds") 

match_numbers=readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_top_motif_matches_human_final_version2.2.Rds")


df=data.frame(ID=names(match_numbers),match_numbers=as.vector(match_numbers))

# Use the merge function to add a column from df2 to df1
merged_representatives <- merge(representatives, df[,c('ID','match_numbers')], by = 'ID', all.x = TRUE)

#Minimum integer threshold used to call the motifs

floor(min(top_motif_matches_human_final[[1]]$score))
min(top_motif_matches_human_final[[1]]$score)

f1 <- function(x) floor(min(x$score))

thresholds=lapply(top_motif_matches_human_final, f1) #, simplity="array"
thresholds_double=lapply(top_motif_matches_human_final, function(x) min(x$score)) #, simplity="array"
# 
df=data.frame(ID=names(top_motif_matches_human_final),
              threshold=unlist(thresholds), threshold_double=unlist(thresholds_double))
# 
 merged_representatives2 <- merge(merged_representatives, df[,c('ID','threshold', 'threshold_double')], by = 'ID', all.x = TRUE)
# 
 write_tsv(merged_representatives2, "../../PWMs_final_version2.2/metadata_representatives_match_numbers_thresholds.tsv") #3294
# 

saveRDS(merged_representatives2,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/metadata_representatives_match_numbers_thresholds.Rds") 

#Motif number distribution
pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/motif-match-number-distribution-final_version2.2.pdf" ))
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
df=data.frame(match_nro=match_numbers)
ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
  scale_y_cut(breaks=c(0.05), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, binwidth 20000')+ theme_minimal()
dev.off()

representative_motif_matches=top_motif_matches_human_final[rep_motifs]
saveRDS(representative_motif_matches, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final_version2.2.Rds")

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

# Add mean PhyloP threshold for conservation

representatives=read_tsv( "/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives_match_numbers_thresholds.tsv") #3294

representatives$phyloP_threshold=NA

files=dir("/scratch/project_2007567/conservation_thresholds_final_version2.2_correct")
for(file in files){
  #file=files[1]
  
  ID=gsub(".csv", "",file)
  thresholds=as.data.frame(read_delim(paste0("/scratch/project_2007567/conservation_thresholds_final_version2.2_correct/", file), delim=" "))
  print(which(representatives$ID==ID))
  representatives$phyloP_threshold[which(representatives$ID==ID)]=thresholds$conservation_threshold_signif
  
}

write_tsv(representatives, "../../PWMs_final_version2.2/metadata_representatives_match_numbers_thresholds_MeanPhyloP_threshold.tsv") #3294
# 

saveRDS(representatives,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/metadata_representatives_match_numbers_thresholds_MeanPhyloP_threshold.Rds") 

#What about the artificial half sites?

halfsites=read_tsv( "/scratch/project_2006203/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv") #3294
halfsites=halfsites[-1,]


top_motif_matches_halfsites=readRDS(file="/scratch/project_2006203/TFBS//ATAC-seq-peaks/RData/top_motif_matches_human_artificialHTSelex_version2.2.Rds") # 
top_motif_matches_halfsites=top_motif_matches_halfsites[-1]

print(length(top_motif_matches_halfsites))

print(head(names(top_motif_matches_halfsites)[which(!(names(top_motif_matches_halfsites) %in% halfsites$ID))]))



match_numbers=unlist(lapply(top_motif_matches_halfsites, length))

df=data.frame(ID=names(match_numbers),match_numbers=as.vector(match_numbers))

# Use the merge function to add a column from df2 to df1
merged_halfsites <- merge(halfsites, df[,c('ID','match_numbers')], by = 'ID', all.x = TRUE)

#Minimum integer threshold used to call the motifs

floor(min(top_motif_matches_halfsites[[1]]$score))
min(top_motif_matches_halfsites[[1]]$score)

f1 <- function(x) floor(min(x$score))

thresholds=lapply(top_motif_matches_halfsites, f1) #, simplity="array"
thresholds_double=lapply(top_motif_matches_halfsites, function(x) min(x$score)) #, simplity="array"
# 
df=data.frame(ID=names(top_motif_matches_halfsites),
              threshold=unlist(thresholds), threshold_double=unlist(thresholds_double))
# 
halfsites <- merge(merged_halfsites, df[,c('ID','threshold', 'threshold_double')], by = 'ID', all.x = TRUE)



halfsites$phyloP_threshold=NA

files=dir("/scratch/project_2007567/conservation_thresholds_final_version2.2_correct")
for(file in paste0(halfsites$ID, ".csv")){
  #file=files[1]
  
  ID=gsub(".csv", "",file)
  thresholds=as.data.frame(read_delim(paste0("/scratch/project_2007567/conservation_thresholds_final_version2.2_correct/", file), delim=" "))
  print(which(representatives$ID==ID))
  halfsites$phyloP_threshold[which(halfsites$ID==ID)]=thresholds$conservation_threshold_signif
  
}

common_columns=intersect(names(representatives), names(halfsites))

test<- rbind(representatives[, common_columns], halfsites[,common_columns])

test$study[test$study=="fromYimeng"]="Current study"
test$study[test$study=="artificial_halfsites"]="Artificial halfsite"

custom_levels <- c("Current study","Jolma2013", "Jolma2015",   "Nitta2015",  "Morgunova2015", "Yin2017", "Artificial halfsite")

# Convert the 'group' column to a factor with the specified levels
test$study <- factor(test$study, levels = custom_levels)

#sort based on study, artificial halfsites last 
test <- test %>%
  arrange(study,experiment)



#"ID"	"family"	"study"	"experiment"	"motif match numbers in GRCh38"	"mean phyloP threshold"

test=test[,c("ID",                
  "Lambert2018_families",
  "study",               
"experiment",
"match_numbers",
"phyloP_threshold") ]

names(test)=c("ID",	"family",	"study",	"experiment",	"motif match numbers in GRCh38",	"mean phyloP threshold")


write_tsv(test, "../../PWMs_final_version2.2/metadata_with_halfsites.tsv") #3294
# 

saveRDS(test,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/metadata_with_halfsites.Rds") 
