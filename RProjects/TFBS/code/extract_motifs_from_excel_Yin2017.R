

#############Read excel tables################

library(tidyverse)
#library(readr)
library(dplyr)
library(tidyverse)
library(readxl)
#BiocManager::install("universalmotif")
library(universalmotif)
library("TFBSTools")
library("Biostrings")
#BiocManager::install("TFBSTools")

rm(list=ls())

source("code/read_excel_tables.R")

# This table contains background subtracted PWMs for all of the TF pairs and individual TFs analyzed.  First line for each factor contains information as follows:	
# Base:	A, C, G or T
# symbol(s):	HGNC symbol (gene names)
# clone type:	Type of construct used: extended DNA binding domain (eDBD) or full-length (FL)
# families:	TF structural family
# experiment type:	Type of experiment: HT-SELEX or methyl-HT-SELEX
# ligand sequence:	Name of primer used to generate selection ligand; number before N indicates length of random sequence
# batch:	batch of HT-SELEX and Methyl-HT-SELEX experiment
# seed:	Sequence used to identify subsequences to be included in the model
# multinomial:	Hamming distance to seed used to identify subsequences to be included in the model
# cycle:	SELEX cycle used to construct the model; previous cycle was used as background unless indicated otherwise by a letter "b" followed by the used background cycle, "u" letter indicates that only the first instance of identical sequencing reads were used in the construction of the PWM
# Matrix is one of the representative PWMs:	value yes or no
# comment:	comment field


all_table <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=22, n_max=8996-22, col_types=c(rep("text",12), rep("numeric",19)))

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,-1]))



colnames(PWMs_metadata)=c("symbol",	"clone","family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment")
PWMs_metadata$organism="Homo_sapiens"
PWMs_metadata$study="Yin2017"
PWMs_metadata$short=NA
PWMs_metadata$type=NA


#symbol     clone family organism     study      experiment ligand        batch  seed               multinomial cycle representative short type              comment filename         

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment")
)

#"ID" "Lambert2018.families"  "new_representative"   "short"                "type"                 "comment"              "filename"             "IC"                  
#"length"  "IC", "consensus"

#HumanTFs-ccbr-TableS1.csv


#remove / from family names
PWMs_metadata=PWMs_metadata %>% 
  mutate(family = gsub("/|\\.\\d+[A-Za-z]+", "_aka_", family))

PWMs_metadata %>%
  count(symbol) 

PWMs_metadata$filename=NA
PWMs_metadata$ID=NA
PWMs_metadata$IC=NA
PWMs_metadata$IC_universal=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=NA


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) ) #1794

#dir.create(paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens"), recursive=TRUE)
#dir.create(paste0("../../PWMs/Yin2017/pwms_space/","Homo_sapiens"), recursive=TRUE)
#dir.create(paste0("../../PWMs/Yin2017/transfac/","Homo_sapiens"), recursive=TRUE)

pwms=paste0("../../PWMs_final/Yin2017/pwms/","Homo_sapiens/")
pwms_space=paste0("../../PWMs_final/Yin2017/pwms_space/","Homo_sapiens/")
transfac_path=paste0("../../PWMs_final/Yin2017/transfac/","Homo_sapiens/")
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(transfac_path, recursive=TRUE)

append=FALSE



#Remove these from metadata
remove=c()

for(m in 1:length(PWMs_list)){
  
  #Empty matrix
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)    
  #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    
  }else{
  
    
  # Write .pfm tab-separated
  filename=paste0(pwms, 
                  paste0(PWMs_metadata[m,
                                       -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                            "representative", "short", "type", "filename","ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                  ".pfm")
  
  ID=paste0(PWMs_metadata[m,
                       -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                            "ID", "IC", "IC_universal","length", "consensus" ))], collapse="_")
  
  
  #Filename and ID should be unique
  if(filename %in% PWMs_metadata$filename){
    print("Warning! Non-unique filenames and IDs")
  }
  
  
  PWMs_metadata$filename[m]=filename
  PWMs_metadata$ID[m]=ID
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=filename, sep="\t")
  
  pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
  colnames(pcm)=NULL
  rownames(pcm)=c("A", "C", "G", "T")
  
  pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                             profileMatrix=pcm
  )
  
  #pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.01)
  icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
  PWMs_metadata[m, "IC"]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
  
  #motif length
  PWMs_metadata[m, "length"]=length(pfm)
  
  motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
  motif@name=paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "representative","short", "type", "filename",
                                                                         "ID", "IC", "IC_universal","length", "consensus"))], collapse="_")
  
  PWMs_metadata[m, "IC_universal"]=motif@icscore
  PWMs_metadata[m, "consensus"]=motif@consensus
  
  #Write transfac format
  transfac=paste0(transfac_path,
                  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone","family","comment", "study","organism","representative","short", "type","filename",
                                                                              "ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),".pfm")
  write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
  #Write .pfm space separated
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pwms_space,
                          paste0(PWMs_metadata[m,
                                               -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative","short", "type", "filename",
                                                                                    "ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                          ".pfm"), sep=" ")
  
  #Write .scpd format, all into a single file
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "representative", "short", "type", "filename",
                                                                                       "ID", "IC", "IC_universal","length", "consensus"))], collapse="_")),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Yin2017/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Yin2017/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  
  
  
  }
}

PWMs_metadata=PWMs_metadata[-remove,]

#Add Family Info



#HumanTFs-ccbr-TableS1.csv
TF_family_mappings <- read_csv("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)
PWMs_metadata$DBD=NULL

table(PWMs_metadata$experiment)


#PWMs_metadata=readRDS(file="RData/Yin2017.Rds")

Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") #923

Methyl_SELEX_motif_categories=read_delim("tables/Methyl_SELEX_metadata_with_info.csv", 
                                                                            delim = ";", 
                                         escape_double = FALSE, trim_ws = TRUE)

PWMs_metadata=PWMs_metadata %>%
  left_join(select(Methyl_SELEX_motif_categories, "ID", "final suggestion" ), by = 'ID')

colnames(PWMs_metadata)[23]="Methyl.SELEX.Motif.Category"

table(PWMs_metadata$Methyl.SELEX.Motif.Category) #some typos

PWMs_metadata$Methyl.SELEX.Motif.Category[PWMs_metadata$Methyl.SELEX.Motif.Category=="No GpG"]="No CpG"
PWMs_metadata$Methyl.SELEX.Motif.Category[PWMs_metadata$Methyl.SELEX.Motif.Category=="Multiple Effects"]="Multiple effects"

categories=table(PWMs_metadata$Methyl.SELEX.Motif.Category)

pie(categories, 
    main = "Proportions of different Methyl-SELEX motif categories", 
    col = rainbow(length(categories)), 
    labels = paste(names(categories), "\n", categories))

remove=c("Little effect", "MethylMinus", "No CpG")

PWMs_metadata=PWMs_metadata[-which(PWMs_metadata$Methyl.SELEX.Motif.Category %in% remove),]

sum(table(PWMs_metadata$Methyl.SELEX.Motif.Category)) #297

#HT-SELEX Methyl-HT-SELEX 
#864             297 

write.table(PWMs_metadata, file="../../PWMs_final/Yin2017/metadata.csv", row.names = FALSE,sep="\t")
saveRDS(PWMs_metadata, file="RData/Yin2017.Rds")


#The categories MethylMinus, MethylPlus, No CpG, Little effect, Multiple effects, Inconclusive were collected using the scripts below
#and manual comparisons to the Yin et al 2017 Science aaj2239_yin_data_s2.pdf Database S2

#From Methyl-selex motifs remove No CpG, MethylMinus and Little effect motifs

#display_table_shape(all_table)



#Sahu et al 2022

#Known TF motifs, a set of X HT-SELEX motifs were collected (refs 10,29 and 30 and unpublished draft motifs). 
#A more compact set of motifs representing different binding specificities was generated by first constructing a dominating set (880)
#covering motifs from the above sources suing the same method and motif distance threshold as in 29.

#In order to remain information of TF binding differences between methylated and non-methylated DNA ligand (10), 
#for each methyl (or non-methyl) motif in the dominating set, the closest non-methyl (methyl) motif of the same TF was added to the
#representative set (respectively) if it was not yet in the set. This resulted in a representative set of 1121 HT-SELEX motifs (Supplementary Table 3).

#In cases where HT-SELEX motifs for several TFs are highly similar, the motifs have been named in figures according to TF class or subclass,
#in order to highlight that we do not know which of the TFs that have similar motifs binds to the motif in the cells. 
#Same principle has been applied also when specificities of closely related TFs have not been measured.
#Supplementary Table 5 shows the naming for the motifs in each figure.

#Additionally, the following promoter core motifs were collected for TSS analyses from literature: 
#TATA box, Initiator, CCAAT-box, GC-box from ref. 31, BRE, MTE, DPE from ref. 32, and BANP from ref. 33.

#A control set of reversed but not complemented motifs was generated by reversing the column order of each motif matrix


#Methyl minus/plus etc information
#Table S3 - bisulfite-SELEX dinucleotide counts and enrichment	

#"symbol",	"clone","family", "organism",	"study","experiment",
#"ligand",	"batch", "seed",	"multinomial",
#"cycle","representative", "short", "type","comment

#section a: bisulfite SELEX analysis of mCG enrichment	
#TF_name	#HGNC name of TF
#clone Type of clone construct used: extended DNA binding domain (eDBD) or full-length (FL)
#call	#final classification of CG positions 
#to Little effect, 
#or MethylPlus 
#or MethylMinus; 
#cases where data is inconsistent are labeled as 'inconclusive', 
#asterisk indicates TFs with different bisulfite-SELEX calls in different experiments. 
#final call is based on bisulfite call except for cases where bisulfite experiment had low enrichment, low complexity seed, or no data.
#methyl-SELEX_call	call from methyl-SELEX
#bisulfite_call	#call from bisulfite-SELEX. Note also that cut-off was set at 10% and 50% to indicate weak effect and strong effect.
#comment	#comment field
#seed	#seed used to count the dinucleotides (counted dinucleotide position is in bold)
#start_position_of_CG	#position of the counted dinucleotide
#enrichment_of_mCG	#enrichment of methyl-CG (see Methods for details)
#CG_frequency	#frequency of CG at the position indicated
#TG_frequency	#frequency of TG at the position indicated



bisulfite_selex <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                              sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=21, n_max=543-21,col_types=c(rep("text",7), rep("numeric",10))) #n_max=8996-22, 

#b table, replicate experiments	
bisulfite_selex_replicate <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                        sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=548, n_max=599-548,col_types=c(rep("text",7), rep("numeric",10)))


colnames(bisulfite_selex)=colnames(bisulfite_selex_replicate)=c("symbol",
                                                                "clone",
                                                                "call",
                                                                "methyl-SELEX_call",
                                                                "bisulfite_call", 
                                                                "comment",
                                                                "seed",
                                                                "start_position_of_CG",
                                                                "enrichment_of_mCG",
                                                                "CG_frequency_cycle_3_normal_seq",
                                                                "TG_frequency_cycle_3_normal_seq",
                                                                "CG_frequency_cycle_3_bisulfite_seq",
                                                                "TG_frequency_cycle_3_bisulfite_seq",
                                                                "CG_frequency_cycle_4_normal_seq",
                                                                "TG_frequency_cycle_4_normal_seq",
                                                                "CG_frequency_cycle_4_bisulfite_seq",
                                                                "TG_frequency_cycle_4_bisulfite_seq")

naind=which(is.na(bisulfite_selex$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  bisulfite_selex[g,"symbol"]=bisulfite_selex[g[1]-1,"symbol"]
  
}

naind=which(is.na(bisulfite_selex_replicate$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  bisulfite_selex_replicate[g,"symbol"]=bisulfite_selex_replicate[g[1]-1,"symbol"]
  
}

#which symbol has asterisk in its name

asterisk_ind=grep("\\*", bisulfite_selex$symbol)
bisulfite_selex$asterisk="NO"
bisulfite_selex$asterisk[asterisk_ind]="YES"
bisulfite_selex$symbol=gsub("\\*","", bisulfite_selex$symbol)

asterisk_ind=grep("\\*", bisulfite_selex_replicate$symbol)
bisulfite_selex_replicate$asterisk="NO"
bisulfite_selex_replicate$asterisk[asterisk_ind]="YES"
bisulfite_selex_replicate$symbol=gsub("\\*","", bisulfite_selex_replicate$symbol)

bisulfite_selex$replicate="NO"
bisulfite_selex_replicate$replicate="YES"

bisulfite_selex_all=rbind(bisulfite_selex, bisulfite_selex_replicate)


#section c: TFs for which bisulfite data was not obtained
#TF_name		call	methyl-SELEX_call		comment	seed for methyl-SELEX

TFs_without_bisulfite_data <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                         sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=605, n_max=703-605,col_types=c(rep("text",7)))

TFs_without_bisulfite_data =TFs_without_bisulfite_data[,-c(2,5)]
colnames(TFs_without_bisulfite_data) <- c("symbol", "call", "methyl-SELEX_call", "comment", "seed") #seed_methyl-SELEX

naind=which(is.na(TFs_without_bisulfite_data$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  TFs_without_bisulfite_data[g,"symbol"]=TFs_without_bisulfite_data[g[1]-1,"symbol"]
  
}





Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX")
#Match seeds in TFs_without_bisulfite_data and Methyl_SELEX_metadata
Methyl_SELEX_metadata <- Methyl_SELEX_metadata %>%
  left_join(select(TFs_without_bisulfite_data, "symbol", "call","methyl-SELEX_call", "comment", "seed" ), by = c('symbol','seed'))


#All possible outcomes
#Table S5 - Contingency table describing the combined effect of the methyl-preference of a TF, 
#methylation of its site, and biological activity of the TF on the methylated state of the bound CpG

#TFs preference: The binding of TFs to DNA is inhibited (MethylMinus), enhanced (MethylPlus) or not affected (Little effect) by CpG methylation.
#state : The state of the cytosine, non-methylated cytosine (C) or methylated cytosine (mC)
#binding: TFs could (YES) or couldn't (No) bind to subsequences with methylated or non-methylated CpG site(s).
#effect on C:  The effect on cytosine after TFs binding to DNA, 
#inducing methylation on normal cytosine (methylation), 
#demethylation of methylated cytosine (demethylation) 
#or keeping the original state (no effect).

possible_outcomes <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                sheet="S5 All possible outcomes", col_names=FALSE, skip=8, n_max=32-8,col_types=c(rep("text",5)))
colnames(possible_outcomes)=c("TFs_preference", "state", "binding", "effect_on_C", "outcome")

#Is the methyl plus/minus info unique for each TF
TFs_unique=unique(bisulfite_selex_all$symbol) #289

methyl_info<-list()

#All methylation info for each TF protein
for(TF in TFs_unique){
  #TF=TFs_unique[1]
  methyl_info[[TF]]=unique(bisulfite_selex_all$call[which(bisulfite_selex_all$symbol==TF)])
  
}

max_length <- max(sapply(methyl_info, length))

# Pad the shorter vectors with NA values
lst_padded <- lapply(methyl_info, function(v) {
  length(v) <- max_length
  return(v)
})

# Combine the vectors by row
result <- do.call(rbind, lst_padded)

#which motifs have unique methylation status or na
na_counts <- apply(result, 1, function(row) sum(is.na(row)))

which(na_counts==3)

df_methyl_info=data.frame(call=result[which(na_counts==3),1])
df_methyl_info$symbol=rownames(df_methyl_info)

Methyl_SELEX_metadata <- Methyl_SELEX_metadata %>%
  left_join(select(df_methyl_info, "symbol", "call", ), by = c('symbol'))


which(na_counts==2)
#KLF17    RFX5    BCL6   BCL6B CREB3L1    E2F2    HSF5     YY1 
#64     114     155     156     162     167     200     232 

which(na_counts==1)
#PAX3 PAX7 
#94   96 
which(na_counts==4)
#empty

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="KLF17"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="KLF17",c("symbol", "clone", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="KLF17",c("call.y")]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="KLF17")[2:1],"call"] 

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="RFX5"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="RFX5",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="RFX5", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="RFX5")[c(1,3)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6")[c(1,2)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6B"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6B",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6B", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6B")[c(1,2)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="CREB3L1",c("symbol", "clone", "seed", "call.y")]

library(DECIPHER)
library(Biostrings)

collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="CREB3L1",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}
for(i in 1:length(methyl_selex)){
  Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol=="CREB3L1")[i], "call.y"]=unique(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),"call"][collect_match_ind[[i]],"call"])
}

# getOption("digits")
# options(digits = 22)
# 
# pairwiseAlignment(bisulfite, methyl_selex[2])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# 
# 
# query = DNAStringSet( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)
# target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )
# names(target)=paste0(Methyl_SELEX_metadata$symbol,"_", Methyl_SELEX_metadata$clone)
# 
# Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol == concordant_TF),]
# 
# alignment=pairwiseAlignment(target, query[1])
# scores=alignment@score
# names(scores)=Methyl_SELEX_metadata$symbol
# 
# 
# scores_order=order(scores, decreasing=TRUE)
# scores_sorted=sort(scores, decreasing=TRUE)


bisulfite_selex_all[which(bisulfite_selex_all$symbol=="E2F2"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="E2F2",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="E2F2", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="E2F2")[c(5,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="HSF5"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="HSF5",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="HSF5", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="HSF5")[c(2,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="YY1"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="YY1",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="YY1", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="YY1")[c(2,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX3"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3",c("symbol", "clone", "seed", "call.y")]

collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX3"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3", "call.y"]=c("Multiple effects", "MethylPlus", "MethylMinus")


bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX7"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7",c("symbol", "clone", "seed", "call.y")]

collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX7"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7", "call.y"]=c("MethylMinus", "Multiple Effects",  
"MethylPlus",
"MethylMinus",  
"Multiple Effects",  
"MethylPlus")

#Test which motifs do not have CpG

Methyl_SELEX_metadata$no_cpg=""


dinucleotides <- c("CG", "YG", "MG", "SG", "CK", "CR", "CS")
                #   "YK", "YR", "YS", "MK", "MR", "MS", "SK", "SR", "SS")

# Function to check for dinucleotides
check_dinucleotides <- function(sequence, dinucleotides) {
  result=FALSE
  for (dinucleotide in dinucleotides) {
    positions <- gregexpr(pattern = dinucleotide, sequence)[[1]]
    if (positions[1] != -1) {
      #print(paste("Dinucleotide", dinucleotide, "found at positions:", paste(positions, collapse=", ")))
      result=TRUE
      
    }
  }
  result
}

cpg_sites=sapply(as.list(Methyl_SELEX_metadata$consensus), check_dinucleotides, dinucleotides=dinucleotides)

Methyl_SELEX_metadata$no_cpg[which(cpg_sites==FALSE)]="No CpG"



write.table(Methyl_SELEX_metadata,file="tables/Methyl_SELEX_metadata_with_info.tsv",sep="\t", row.names = FALSE)



#Below this some trials, not used
stop()


#Are the TFs in the replicate table subset of the main table
TFs_unique_replicate=unique(bisulfite_selex_replicate$symbol) #30
length(which(TFs_unique_replicate %in% TFs_unique)) #28

TFs_unique_replicate[which(!(TFs_unique_replicate %in% TFs_unique))] #all are included, note the asterisk

#Are the TFs in the TFs_without_bisulfite_data missing from the bisulfite selex table
TFs_unique_without_bisulfite_data=unique(TFs_without_bisulfite_data$symbol) #82
length(which(TFs_unique_without_bisulfite_data %in% TFs_unique)) #0

#match metadata and methylation info
as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1] #541

Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX")
write.table(Methyl_SELEX_metadata, file="tables/Methyl_SELEX_metadata.tsv",sep="\t", row.names = FALSE)
write.table(bisulfite_selex,file="tables/bisulfite_selex.tsv",sep="\t", row.names = FALSE)
write.table(bisulfite_selex_replicate,file="tables/bisulfite_selex_replicate.tsv",sep="\t", row.names = FALSE)
write.table(TFs_without_bisulfite_data,file="tables/TFs_without_bisulfite_data.tsv",sep="\t", row.names = FALSE)

methyl_metadata_names=c("methyl_symbol",
                        "methyl_clone",
                        "methyl_call",                              
                        "methyl-SELEX_call",
                        "bisulfite_call",
                        "methyl_comment",                           
                        "methyl_seed",
                        "start_position_of_CG",
                        "enrichment_of_mCG",                 
                        "CG_frequency_cycle_3_normal_seq",
                        "TG_frequency_cycle_3_normal_seq",  
                        "CG_frequency_cycle_3_bisulfite_seq",
                        "TG_frequency_cycle_3_bisulfite_seq",
                        "CG_frequency_cycle_4_normal_seq",    "TG_frequency_cycle_4_normal_seq",   
                        "CG_frequency_cycle_4_bisulfite_seq", "TG_frequency_cycle_4_bisulfite_seq")

methyl_metadata=as.data.frame(matrix("", nrow(Methyl_SELEX_metadata), length(methyl_metadata_names)))
names(methyl_metadata)=methyl_metadata_names

Methyl_SELEX_metadata=cbind(Methyl_SELEX_metadata, methyl_metadata)

concordant=as_tibble(bisulfite_selex[which(bisulfite_selex$comment=="concordant"),])

bisulfite_selex
target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )

collect_match_ind=c()

Methyl_SELEX_metadata$seed_match_score=NA

#TFs_unique=unique(concordant$symbol)
#TFs_unique=TFs_unique[-which(is.na(TFs_unique))]

TFs_unique=unique(bisulfite_selex$symbol)
TFs_unique=TFs_unique[-which(is.na(TFs_unique))]

for(concordant_TF in TFs_unique){
  #concordant_TF=unique(concordant$symbol)[1]
  print(concordant_TF)
  
  #concordant_table=concordant %>% filter(symbol==concordant_TF)
  #query = DNAStringSet( as.vector(concordant_table %>% dplyr::select(seed))$seed)
  
  concordant_table=bisulfite_selex %>% filter(symbol==concordant_TF)
  
  query = DNAStringSet( as.vector(concordant_table %>% dplyr::select(seed))$seed)
  

  for(i in 1:length(query)){
    print(i)
    q=query[i]
    
    alignment=pairwiseAlignment(target, q)
    scores=alignment@score
    names(scores)=Methyl_SELEX_metadata$symbol
    scores_order=order(scores, decreasing=TRUE)
    scores_sorted=sort(scores, decreasing=TRUE)
    scores_sorted[which(names(scores_sorted)==concordant_TF)]
    
    
    match_ind=scores_order[which(names(scores_sorted)==concordant_TF)][1]
    if(!is.na(match_ind)){
      if(scores_sorted[which(names(scores_sorted)==concordant_TF)][1]>0){ #Score needs to be higher than 0
      
      #Is the match_ind already filled
      if(!(match_ind %in% collect_match_ind)){
      
      #Do the clone match
        if(Methyl_SELEX_metadata$clone[match_ind]==concordant_table$clone[i]){
          Methyl_SELEX_metadata$seed_match_score[match_ind]=scores_sorted[which(names(scores_sorted)==concordant_TF)][1] #score
          Methyl_SELEX_metadata[match_ind, names(methyl_metadata)]=concordant_table[i,]
      
        }else{
          print("clone does not match")
        }
      }else{
        print("already filled")
      }
    
      collect_match_ind=c(collect_match_ind, match_ind)
 
      }
    }
  }
  
}

write.table(Methyl_SELEX_metadata, file="tables/Methyl_SELEX_metadata.tsv",sep="\t", row.names = FALSE)


concordant_TF=unique(concordant$symbol)[1]
seq1 = DNAString(as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)

BiocManager::install("DECIPHER")

library(DECIPHER)

collect_match_ind=c()

sequence_list=c( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed, as.vector(Methyl_SELEX_metadata$seed))

sequences = DNAStringSet( sequence_list )

query = DNAStringSet( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)
target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )
names(target)=paste0(Methyl_SELEX_metadata$symbol,"_", Methyl_SELEX_metadata$clone)

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol == concordant_TF),]

alignment=pairwiseAlignment(target, query[1])
scores=alignment@score
names(scores)=Methyl_SELEX_metadata$symbol


scores_order=order(scores, decreasing=TRUE)
scores_sorted=sort(scores, decreasing=TRUE)

scores_order[which(names(scores_sorted)==concordant_TF)] #168(this is already taken), 169

match_ind=scores_order[which(names(scores_sorted)==concordant_TF)][1]

collect_match_ind=c(collect_match_ind, match_ind)

Methyl_SELEX_metadata[match_ind), 
                      names(methyl_metadata)]=concordant[]


alignment=pairwiseAlignment(target, query[2])
scores=alignment@score
names(scores)=Methyl_SELEX_metadata$symbol


scores_order=order(scores, decreasing=TRUE)
scores_sorted=sort(scores, decreasing=TRUE)

#Remove those that are already in the collect_match_ind

match_ind=scores_order[which(names(scores_sorted)==concordant_TF)] #168(this is already taken), 169

match_ind=match_ind[-which(match_ind %in% collect_match_ind)][1]

collect_match_ind=c(collect_match_ind, match_ind)




methyl_minus=bisulfite_selex[which(bisulfite_selex$comment=="concordant" & bisulfite_selex$call=="MethylMinus"),]

methyl_minus_TF=unique(methyl_minus$symbol)

which(!(Methyl_SELEX_metadata$symbol %in% methyl_minus_TF))

for( TF in methyl_minus_TF){
  TF=methyl_minus_TF[3]
  print(length(which(Methyl_SELEX_metadata$symbol==TF) )) #This can be 0,1,2,3,4
  #If two, they are from different clones
}

methyl_minus[methyl_minus$symbol==TF,c("seed","start_position_of_CG")]
methyl_minus[methyl_minus$symbol==TF,]

NRTGAYGTYAYN                    6
NRTGAYGTGGYN                    6
NGCCACRTCAYN YKCCACRTCAYC
NRTGAYGTCAYN GRTGAYGTCAYN
NGCCACGTGKMN NGCCACGTGKMN


seq1 = DNAString("NRTGAYGTYAYN")
seq2 = DNAString("NRTGAYGTGGYN")

seq3 = DNAString("NGCCACRTCAYN")
seq4 = DNAString("NRTGAYGTCAYN")
seq5 = DNAString("NGCCACGTGKMN")

# Perform pairwise alignment
pairwiseAlignment(seq1, seq3)
pairwiseAlignment(seq1, seq4) #Highest score
pairwiseAlignment(seq1, seq5)

pairwiseAlignment(seq2, seq3)
pairwiseAlignment(seq2, seq4) #Highest score
pairwiseAlignment(seq2, seq5)


# Print alignment
print(alignment)
print(alignment2)


# Calculate percentage identity
matches = sum(matchPattern(alignment, alignment))
total = width(alignment)
percentage_identity = matches / total * 100

print(paste("Percentage identity: ", percentage_identity))



bisulfite_selex %>% filter(symbol==TF,) #FL eDBD
PWMs_metadata[which(PWMs_metadata$symbol==TF),]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),c("seed", "consensus")] #eDBD FL
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] #eDBD FL

TF=methyl_minus_TF[2]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=methyl_minus[methyl_minus$symbol==TF,]

TF=methyl_minus_TF[3] #Both tables have two hits, eDBD and FL, sort methyl_minus based on the order in metadata
Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)]

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), 
                      names(methyl_metadata)]=methyl_minus[which(methyl_minus$symbol==TF)[match(Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)], 
                                                                                                methyl_minus$clone[methyl_minus$symbol==TF] )],]
#Same info twice, the seeds are equal at the interesting positions, is this even CpG
TF=methyl_minus_TF[4]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=rbind(methyl_minus[methyl_minus$symbol==TF,],methyl_minus[methyl_minus$symbol==TF,])

TF=methyl_minus_TF[5] #Same situation as above
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=rbind(methyl_minus[methyl_minus$symbol==TF,],methyl_minus[methyl_minus$symbol==TF,])

TF=methyl_minus_TF[6]
#Which seeds in metadata are of same length?


same_length=Methyl_SELEX_metadata$length[which(Methyl_SELEX_metadata$symbol==TF)]==unique(nchar(methyl_minus[methyl_minus$symbol==TF,"seed"]$seed))

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF)[which(same_length)],names(methyl_metadata)]=methyl_minus[which(methyl_minus$symbol==TF)[match(Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)], 
                                                                                                                                      methyl_minus$clone[methyl_minus$symbol==TF] )],]
TF=methyl_minus_TF[7]


#GTYAYGTGAY

NRTGAYGTYAYN
NRTGAYGTGGYN
PWMs_list[[which(PWMs_metadata$symbol==methyl_minus_TF)[1]]]
PWMs_list[[which(PWMs_metadata$symbol==methyl_minus_TF)[2]]]



PWMs_metadata[which(PWMs_metadata$symbol==methyl_minus_TF),]

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==methyl_minus_TF),]
#RTCAYGTGMN
#NNNNCGNNNN

NGCCACRTCAYN
NRTGAYGTCAYN
NGCCACGTGKMN


TF=unique((PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol))$symbol)[5]

#IUPAC nucleotide codes https://www.bioinformatics.org/sms/iupac.html


Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==bisulfite_selex$symbol[1]),"seed"]
"CYAATTAN"
"NTCGTTAN"
"CYAATTAN"
"NTCGTTAN"

bisulfite_selex[1,"seed"]
TTAAYGAW

bisulfite_selex[3,]
NGTYGTAAAWN
bisulfite_selex[3,"seed"]
NGTYGTAAAAN

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==bisulfite_selex$symbol[2]),]
"NGTCGTAAANN"
"NGYAATAAAAN"


for(TF in unique((PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol))$symbol)[1:5]){
  print("TF")
  print(TF)
  print("Methyl_SELEX")
  print(Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),])
  print(Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"])
  print("bisulfite_SELEX")
  print(bisulfite_selex[which(bisulfite_selex$symbol==TF),])
  print(bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"])
  print("bisulfite_SELEX_replicate")
  print(bisulfite_selex_replicate[which(bisulfite_selex_replicate$symbol==TF),])
  print(bisulfite_selex_replicate[which(bisulfite_selex_replicate$symbol==TF),"seed"])
  print("TFs_without_bisulfite_data")
  print(TFs_without_bisulfite_data[which(TFs_without_bisulfite_data$symbol==TF),])
  #  print(TFs_without_bisulfite_data[which(TFs_without_bisulfite_data$symbol==TF),"seed"])
  
  #PROX1: seed matches with bisulfite selex
  
  tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] 
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed 
  if(seed1==seed2) #Seed does not match
    
    #MSX2
    tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] #two
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed #only one
  if(seed1==seed2) #Seed does not match
    
    #MSX1
    tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] #two
  tmp2=bisulfite_selex[which(bisulfite_selex$symbol==TF),]
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed #only one
  if(seed1==seed2) #Seed does not match
    
}

PWMs_metadata[which(PWMs_metadata$symbol==as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1]), ]
bisulfite_selex[which(bisulfite_selex$symbol==as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1]), "seed"]

test=PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% left_join(bisulfite_selex, by=c("seed"))

flights2 %>% left_join(planes, by = "tailnum")


which(PWMs_metadata$experiment=="Methyl-HT-SELEX")


