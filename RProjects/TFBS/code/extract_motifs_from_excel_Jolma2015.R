

#############Read excel tables################

library(tidyverse)
library(readxl)
#all_table <- read_excel("useful/tidyverse_notes/utility/multiple_tables_sheet.xlsx", col_names=FALSE)



#Jolma2015
rm(list=ls())


source("code/read_excel_tables.R")
# This table contains background subtracted PWMs for all of the TF pairs and individual TFs analyzed.  
# First line for each factor contains information as follows:	
# Base:	A, C, G or T
# symbol(s):	HGNC symbol(s) (gene names); Either for reference models that were generated for individual 
  # TFs in this study or for the co-operative pairs, in which case the TF names are separated with an underscore "_". 
  # First of the TFs is TF1 that was expressed as a SBP tagged clone and was used in first of the affinity separations 
  # and the second one is the 3xFLAG tagged TF2
# families:	TF structural families for the TFs. Families of the two TFs in a pair are separated by an underscore
# experiment type:	Type of experiment: HT-SELEX or CAP-SELEX
# ligand sequence:	Name of primer used to generate selection ligand; number before N indicates length of random sequence
# batch:	batch of HT-SELEX or CAP-SELEX experiment
# seed:	Sequence used to identify subsequences to be included in the model
# multinomial:	Hamming distance to seed used to identify subsequences to be included in the model
# cycle: 	SELEX cycle used to construct the model; previous cycle was used as background unless indicated otherwise by a letter 
  # "b" followed by the used background cycle, "u" letter indicates that only the first instance of identical sequencing reads 
  # were used in the construction of the PWM
# comment: 	comment field
# Matrix is one of the representative PWMs	value: yes or no
# Genomic conservation adjusted p-value:	Adjusted p-value for the conservation of the PWM target sites in comparison to target sites 
  # located on control regions. See methods for details
# Enrichment in LoVo ChIP-seq clusters p-value:	p-value for the enrichment of CAP-SELEX pwm target sites over control regions. See methods for details
# Category:	Whether the model is both enriched in the ChIP-seq clusters (See Yan et al 2013) and evolutionarily conserved in the 
  # genome "enriched and conserved", either "enriched or conserved", or "neither"
# Models marked in orange (last dimer models) are technical replicates performed either with CAP-SELEX (DBD+ clones) or with HT-SELEX (Full length TFs), 
# models marked in lila are individual HT-SELEX models described for the first time in this study 
# and models marked with brown are enriched motifs from ChIP-Exo experiments	

oranges=which(seq(20, 3335,5) %in% seq(2830,3160,5)) #67
lilas=which(seq(20, 3335,5) %in% seq(3165,3325,5)) #33
browns=which(seq(20, 3335,5) %in% seq(3330,3335,5)) #2

all_table <- read_excel("../../PWMs/Jolma2015/41586_2015_BFnature15518_MOESM33_ESM.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=19, col_types=c(rep("text",13), rep("numeric",20)))
#

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,-1]))

#remove ChIP-EXO motifs
PWMs_metadata=PWMs_metadata[-browns,]


colnames(PWMs_metadata)=c("symbol",	"family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment",
                          "conservation p-value", "ChIP-seq cluster p-value", "category")


#"symbol",	"clone","family", "organism",	"study","experiment",
#"ligand",	"batch", "seed",	"multinomial",
#"cycle","representative", "short", "type","comment", "filename"

PWMs_metadata$clone=NA
PWMs_metadata$study="Jolma2015"
PWMs_metadata$short=NA
PWMs_metadata$type=NA
PWMs_metadata$filename=NA
PWMs_metadata$organism="Homo_sapiens"

PWMs_metadata$representative=toupper(PWMs_metadata$representative)

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
  "ligand",	"batch", "seed",	"multinomial",
  "cycle","representative", "short", "type","comment", "filename"))

#add orange and lila information
na_ind=which(is.na(PWMs_metadata$comment[oranges]))

PWMs_metadata$comment[oranges[na_ind]]="technical-replicates"
PWMs_metadata$comment[oranges[-na_ind]]=paste0(PWMs_metadata$comment[oranges[-na_ind]], ", technical-replicates")

PWMs_metadata$comment[lilas]=paste0(PWMs_metadata$comment[lilas], ", new-individual-models")


#no representative info
ind=which(!PWMs_metadata$representative %in% c("YES", "NO")) # lilas and ind are the same
PWMs_metadata$representative[ind]=NA

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

#remove ChIP-EXO motifs
PWMs_list=PWMs_list[-browns]

dir.create(paste0("../../PWMs/Jolma2015/pwms/","Homo_sapiens"), recursive = TRUE)

append=FALSE

for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Jolma2015/pwms/",PWMs_metadata[m,"organism"],"/", 
                          paste0(PWMs_metadata[m,
                                               -which(colnames(PWMs_metadata)%in% c("clone","family","comment", "study","organism","short", "type","filename"))], collapse="_"),
                          ".pfm"), sep="\t")
  
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], collapse="_")),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Jolma2015/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Jolma2015/",PWMs_metadata[m,"organism"],"_all", ".scpd"))
  
  
  PWMs_metadata$filename[m]=paste0("PWMs/Jolma2015/pwms/Homo_sapiens/", paste0(PWMs_metadata[m,
                                                                                           -which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "short", "type", "filename"))], collapse="_"),".pfm")
  
}

#there are some PWMs twice in the metadata table, the PWM for the latter occurrence is used
PWMs_metadata[which( duplicated(PWMs_metadata$filename)==TRUE),"filename"]

PWMs_metadata[ which(PWMs_metadata$filename %in% as.vector(PWMs_metadata[which( duplicated(PWMs_metadata$filename)==TRUE),]$filename)), ]

write.table(PWMs_metadata, file="../../PWMs/Jolma2015/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(PWMs_metadata, file="data/Jolma2015.Rds")

stop()

PWMs_metadata %>%
  count(symbol) 
#435

PWMs_metadata %>%
  count(family) 


#Number of pairs
PWMs_metadata %>%
  count(symbol)%>%
  filter(str_detect(symbol, "_") )
#318   27 single



#filter oranges and greens
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  count(symbol)
#463



#How many types unique(PWMs_metadata$`site type`)
#[1] "monomeric"             "dimeric"               "monomer"               "monomeric or dimeric"  "dimer of dimers"       "putative multimer"    
#[7] "monomer or dimer"      "trimeric"              "putatively multimeric"

#Number of human TFs representative
PWMs_metadata %>%
  filter(`Matrix is one of the representative PWMs`=="no")%>%
  count(symbol)

#representative 223
#nonrepresentative 167







#display_table_shape(all_table)

