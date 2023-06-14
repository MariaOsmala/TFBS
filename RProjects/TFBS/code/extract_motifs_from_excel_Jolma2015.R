

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
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







colnames(PWMs_metadata)=c("symbol",	"family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment",
                          "conservation p-value", "ChIP-seq cluster p-value", "category")

#Remove Hoxa10 Pantheropsis Guttatus, snake
snake_ind=which(PWMs_metadata$symbol=="Hoxa10")


#remove ChIP-EXO motifs
PWMs_metadata=PWMs_metadata[-c(snake_ind,oranges,browns),]

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
#na_ind=which(is.na(PWMs_metadata$comment[oranges]))
#PWMs_metadata$comment[oranges[na_ind]]="technical-replicates"
#PWMs_metadata$comment[oranges[-na_ind]]=paste0(PWMs_metadata$comment[oranges[-na_ind]], ", technical-replicates")
#PWMs_metadata$comment[lilas]=paste0(PWMs_metadata$comment[lilas], ", new-individual-models")

PWMs_metadata$ID=""
PWMs_metadata$IC=NA
PWMs_metadata$IC_universal=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=""

#no representative info
ind=which(!PWMs_metadata$representative %in% c("YES", "NO")) # lilas and ind are the same
PWMs_metadata$representative[ind]=NA

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

#remove ChIP-EXO motifs and replicates
PWMs_list=PWMs_list[-c(snake_ind,oranges,browns)]

pwms="../../PWMs_final/Jolma2015/pwms/Homo_sapiens"
pwms_space="../../PWMs_final/Jolma2015/pwms_space/Homo_sapiens"
pwms_transfac="../../PWMs_final/Jolma2015/transfac/Homo_sapiens"

dir.create(pwms, recursive = TRUE)
dir.create(pwms_space, recursive = TRUE)
dir.create(pwms_transfac, recursive = TRUE)

append=FALSE

#[1] "Warning! Non-unique filenames and IDs"
#[1] "POU2F1_ETV4_CAP-SELEX_TAGAAC40NAAT_AX_NCCGGATATGCAN_1_2"
#[1] "Warning! Non-unique filenames and IDs"
#[1] "POU2F1_SOX15_CAP-SELEX_TCGCCA40NCCT_AS_WTGMATAACAATR_1_2"
#[1] "Warning! Non-unique filenames and IDs"
#[1] "ERF_DLX2_CAP-SELEX_TCGAAA40NAAT_AAC_RSCGGAANNNNNYMATTA_1_3"
#[1] "Warning! Non-unique filenames and IDs"
#[1] "POU2F1_GSC2_CAP-SELEX_TGCCAT40NACC_AS_NNGATTANNNATGCAWNNN_1_2"

remove=c()

for(m in 1:length(PWMs_list)){
  
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
  }else{
  
    # Write .pfm tab-separated
    filename=paste0(pwms,"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                              "representative", "short", "type", "filename","ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                 "ID", "IC", "IC_universal","length", "consensus" ))], collapse="_")
    
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
      print(ID)
      filename=paste0(pwms,"/", 
                      paste0(PWMs_metadata[m,
                                           -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                                "representative", "short", "type", "filename","ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                      "_v2",".pfm")
      ID=paste0( paste0(PWMs_metadata[m,
                              -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                   "ID", "IC", "IC_universal","length", "consensus" ))], collapse="_"), "_v2")
      
      
      
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,file=filename, sep="\t")
  
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
    motif@name=ID
    
    PWMs_metadata[m, "IC_universal"]=motif@icscore
    PWMs_metadata[m, "consensus"]=motif@consensus
    
    transfac=paste0(pwms_transfac,"/", ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE, 
                file=paste0(pwms_space,"/", ID, ".pfm"), sep=" ")
  
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Jolma2015/all", ".scpd"))
    append=TRUE
  
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Jolma2015/","all", ".scpd"))
  
  
  }
}





#there are some PWMs twice in the metadata table, the PWM for the latter occurrence is used, NOT ANYMORE
#PWMs_metadata[which( duplicated(PWMs_metadata$filename)==TRUE),"filename"]

#PWMs_metadata[ which(PWMs_metadata$filename %in% as.vector(PWMs_metadata[which( duplicated(PWMs_metadata$filename)==TRUE),]$filename)), ]

#Add family mappings

TF_family_mappings <- read_delim("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv", delim=",")

#commondf <- merge(datas, TF_family_mappings, by = 'symbol', all.x = TRUE)

#There is no singles in this data
single <- PWMs_metadata$symbol[str_detect(PWMs_metadata$symbol, "_") == FALSE] #31
single_index <- which(PWMs_metadata$symbol %in% single)
Family_names <- rep(NA, nrow(PWMs_metadata))
Family_names[single_index]=as.vector(PWMs_metadata[single_index,"symbol"] %>%
  left_join(select(TF_family_mappings, "symbol", "DBD" ), by = 'symbol') %>% select("DBD"))$DBD


pair <- PWMs_metadata$symbol[str_detect(PWMs_metadata$symbol, "_") == TRUE]
pair_index <- which(PWMs_metadata$symbol %in% pair) #562
#[pair_index]
pairs <- as.data.frame(do.call(rbind, strsplit(PWMs_metadata$symbol[pair_index], split="_")))
colnames(pairs)=c("symbol", "symbol1")
pairs$order=seq(nrow(pairs))

commondf_pairs <- merge(pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,5)]
commondf_pairs <- commondf_pairs[order(commondf_pairs$order), ]
commondf_pairs$order <- NULL

names(commondf_pairs)=c("symbol0", "symbol", "DBD0")
commondf_pairs$order=seq(nrow(commondf_pairs))
commondf_pairs2 = merge(commondf_pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,4,6)]
commondf_pairs2 <- commondf_pairs2[order(commondf_pairs2$order), ]
commondf_pairs2$order <- NULL
commondf_pairs2 <- commondf_pairs2[, c(2,1,3,4)]
names(commondf_pairs2)=c("symbol0", "symbol1", "DBD0", "DBD1")

pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

Family_names[pair_index]=pair_families


PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = Family_names, .before = 4)
unique_Lambertfamilies <- unique(PWMs_metadata$Lambert2018_families)

#Is some of the Lamber2018_families missing

which(is.na(PWMs_metadata$Lambert2018_families))

PWMs_metadata[582, 1:5]
PWMs_metadata$Lambert2018_families[582]=PWMs_metadata$family[582]

write.table(PWMs_metadata, file="../../PWMs_final/Jolma2015/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(PWMs_metadata, file="Rdata/Jolma2015.Rds")

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

