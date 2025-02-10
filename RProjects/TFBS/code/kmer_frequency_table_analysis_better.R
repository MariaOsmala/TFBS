library("tidyverse")
library(readr)



frequency_table <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/6mer_counts_fromDrive_28082024.csv", 
                                                                          delim = ";", escape_double = FALSE, trim_ws = TRUE) #3914 rows

frequency_table=frequency_table %>% unite(ID, c(TF1_TF2, Barcode, Batch),sep="_", remove=FALSE) 

#frequency_table$TF1_TF2[grep("NKX",frequency_table$TF1_TF2)]

#There are 1394 individual tables for 1336 number of distinct pairs.
length(unique(frequency_table$TF1_TF2)) #1336

#What if the order does not matter, still 1336
#tmp=as.data.frame(do.call(rbind, strsplit(frequency_table$TF1_TF2, "_")))[, c(2,1)] 
#length(unique(c(paste0(tmp$V1, "_", tmp$V2), frequency_table$TF1_TF2) )) #TRUE) #1336

length(unique(frequency_table$ID)) #1394

experiments=unique(frequency_table$ID) #1394

# In total 2198 screened TF-TF pairs displayed specific interaction including 
# 1329 spacing and orientation preferences and 1131 composite motifs. 
# Table S3 TF interaction matrix. Extract those TF pairs for which the interaction 
#belongs to class 2: The interaction between TFs results in spacing and/or orientation preference.

interaction_table_S3 <- as.data.frame(read_delim("~/projects/TFBS/PWMs_final_version2.2/interaction_table_S3.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE))


rownames(interaction_table_S3) <- interaction_table_S3[,1]
interaction_table_S3 <- interaction_table_S3[,-1]

unique(unlist(apply(interaction_table_S3,1,unique )))
# "self" "0"    "2"    "1"    "1,2" 

#test that the matrix is symmetrix wrt. the diagonal
all(interaction_table_S3 == t(interaction_table_S3)) #TRUE

# Consider only the upper triangle  
interaction_table_S3[lower.tri(interaction_table_S3)] <- NA

#Count the number of 0 1 2 1,2

table(unlist(interaction_table_S3))

#0        1     1,2     2     self 
#40580   869    33     1296   293 

#869+33=902 (there can be multiple composites for the same pair --> 1131)
#1296+33=1329

#Find rows and columns with value 2
spacing_interactions <- which(interaction_table_S3 == 2, arr.ind = TRUE) #1296+33=1329

spacing_pairs=data.frame(TF1=rownames(interaction_table_S3)[spacing_interactions[,1]],
                 TF2=colnames(interaction_table_S3)[spacing_interactions[,2]])

composite_interactions <- which(interaction_table_S3 == 1, arr.ind = TRUE) #869+33=902
composite_pairs=data.frame(TF1=rownames(interaction_table_S3)[composite_interactions[,1]],
                 TF2=colnames(interaction_table_S3)[composite_interactions[,2]])

both_types_of_interactions<- which(interaction_table_S3 == "1,2", arr.ind = TRUE) #33

spacing_pairs=rbind(spacing_pairs, data.frame(TF1=rownames(interaction_table_S3)[both_types_of_interactions[,1]],
                 TF2=colnames(interaction_table_S3)[both_types_of_interactions[,2]])) #1329

composite_pairs=rbind(composite_pairs, data.frame(TF1=rownames(interaction_table_S3)[both_types_of_interactions[,1]], #902
                                                TF2=colnames(interaction_table_S3)[both_types_of_interactions[,2]]))


length(unique(frequency_table$TF1_TF2)) #1336

table(unique(frequency_table$TF1_TF2) %in% paste0(spacing_pairs$TF1,"_", spacing_pairs$TF2) )#TRUE
#FALSE  TRUE 
#586   750 
table(unique(frequency_table$TF1_TF2) %in% paste0(spacing_pairs$TF2,"_", spacing_pairs$TF1) )#TRUE
#FALSE  TRUE 
#750   586 
#total 1336

#Other way around
ind1=which(paste0(spacing_pairs$TF1,"_", spacing_pairs$TF2) %in% unique(frequency_table$TF1_TF2) )#TRUE
#740
ind2=which(paste0(spacing_pairs$TF2,"_", spacing_pairs$TF1) %in% unique(frequency_table$TF1_TF2) )#TRUE
#586

#sum 740+586)=1326

#Are some indexes the same

length(intersect(ind1, ind2)) #7
length(union(ind1, ind2)) #1329

#There is no monomer for all TFs, remove those pairs that have these
missing_monomers=c("HHEX","BHLHB8",
"PROX2",
"NKX2-4",
"DBX2",
"BHLHB5",
"ZBTB33",
"PBX2",
"FOXI2",
"gntR",
"NKX2-6")

which(spacing_pairs$TF1 %in% missing_monomers) %in% which(spacing_pairs$TF2 %in% missing_monomers) #All false
which(spacing_pairs$TF2 %in% missing_monomers) %in% which(spacing_pairs$TF1 %in% missing_monomers) #All false

spacing_pairs=spacing_pairs[-which((spacing_pairs$TF1 %in% missing_monomers) | (spacing_pairs$TF2 %in% missing_monomers) ),]#FALSE

table(spacing_pairs$TF1 %in% monomers$symbol)
table(spacing_pairs$TF2 %in% monomers$symbol)

#For which interacting pairs, the k-mer frequency table is missing
ind1=which( !(paste0(spacing_pairs$TF1,"_", spacing_pairs$TF2) %in% unique(frequency_table$TF1_TF2) ) & !(paste0(spacing_pairs$TF2,"_", spacing_pairs$TF1) %in% unique(frequency_table$TF1_TF2) ))#TRUE
length(ind1) #0

#Composites
table(unique(frequency_table$TF1_TF2) %in% paste0(composite_pairs$TF1,"_", composite_pairs$TF2) )#TRUE
#FALSE  TRUE 
#1314    22
table(unique(frequency_table$TF1_TF2) %in% paste0(composite_pairs$TF2,"_", composite_pairs$TF1) )#TRUE
#FALSE  TRUE 
#1325    11 

#Other way around
ind1=which(paste0(composite_pairs$TF1,"_", composite_pairs$TF2) %in% unique(frequency_table$TF1_TF2) )#TRUE
#740
ind2=which(paste0(composite_pairs$TF2,"_", composite_pairs$TF1) %in% unique(frequency_table$TF1_TF2) )#TRUE

length(intersect(ind1, ind2)) #0
length(union(ind1, ind2)) #33


source("code/half_site_recognition_functions.R", echo=FALSE)

metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

#Do I have motifs for the spacing pairs and composite pairs

#unique_symbols_spacing=metadata %>% filter(study=="fromYimeng"& experiment=="CAP-SELEX" & type=="spacing") %>% select(symbol) %>% unique() %>% pull(symbol)#192


#Does all the spacing motifs have the corresponding k-mer frequency table, yes!
#tmp=metadata %>% filter(study=="fromYimeng"& experiment=="CAP-SELEX" & type=="spacing") %>% select(symbol, ligand, batch) %>% mutate(comb=paste0(symbol,"_",ligand ,"_",batch) )
#length(unique(tmp$comb)) #193

#which(tmp %>% pull(comb) %in% frequency_table$ID #TRUE)

#How many of the frequency_table pairs have a motif
#ind1=which(unique(frequency_table$TF1_TF2) %in% (metadata %>% filter(study=="fromYimeng"& experiment=="CAP-SELEX" & type=="spacing") %>% select(symbol) %>% unique() %>%pull(symbol)))
#192

#3
#ind2=which(unique(paste0(tmp$TF2, "_", tmp$TF1)) %in% (metadata %>% filter(study=="fromYimeng"& experiment=="CAP-SELEX" & type=="spacing") %>% select(symbol) %>% unique() %>%pull(symbol)))

#Metadata for artificial half sites, 7 motifs
metadata_halfsites <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

metadata_halfsites$new_representative="NO"
metadata_halfsites$filename=paste0("../../",metadata_halfsites$filename)
metadata_all=rbind(metadata[, colnames(metadata_halfsites)], metadata_halfsites)


monomers=metadata_all %>% filter(experiment=="HT-SELEX") # 1745

#BHLHB2 is not dimeric (homodimer)

monomers$type[grep("BHLHB2", monomers$symbol, ignore.case=TRUE)]=""
monomers$type[grep("MEF2A", monomers$symbol, ignore.case=TRUE)]=""
monomers$type[grep("TCF3", monomers$symbol, ignore.case=TRUE)]=""

experiments_table=data.frame(experiment=experiments, number_of_orientations=NA, TF1="", TF2="",
                             kmer1="", motif1="", orientation1="", start1="",end1="", 
                             kmer2="", motif2="", orientation2="", start2="", end2="")


#Seems that in this table, if there are 4 orientations, the order is
#HT
#TH
#TT
#HH

#If there are two orientations, they are
#TH
#HT

#If there is only one orientation, it is
#TH check this!!!

options(warn=1)
i=1
#which(experiments=="POU3F2_MEIS3_TTGTCC40NTAG_YRIIII")
for(ex in experiments){ #1394
  #ex=experiments[i]
  print(ex)
  data=frequency_table %>% filter(ID==ex) %>% select(-ID)
  
  print(paste0( nrow(data), "_orientations"))
  
  experiments_table[i, "number_of_orientations"]=as.numeric(nrow(data))
  
  TF1=unlist(strsplit(data$TF1_TF2[1], "_"))[1]
  TF2=unlist(strsplit(data$TF1_TF2[1], "_"))[2]
  
  
  experiments_table$TF1[i]=TF1
  experiments_table$TF2[i]=TF2
  
  if(nrow(data)==4 ){
    kmer1=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][1] #"TAAACA"
    kmer2=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][2] #"TGCGGG"
    
  }else{
    #nrow(data)==2 | nrow(data)==1
    kmer1=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][2] #"TAAACA"
    kmer2=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][1] #"TGCGGG"
  }
  

  
  experiments_table$kmer1[i]=kmer1
  experiments_table$kmer2[i]=kmer2
  
  #Remove those with type:  "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric","putatively multimeric"
  #Alternative: Choose the motif with best match to the k-mer
  #If still several: Choose the shortest motif
  #If still several: Choose the highest IC motif
  
  #No motif for these: "BHLHB8", "PROX2"
  
  motif1=try(find_monomer_match_kmer(monomer=TF1, monomers, six_mer=kmer1), silent=TRUE)
  #print(motif1)
  motif2=try(find_monomer_match_kmer(monomer=TF2, monomers, six_mer=kmer2), silent=TRUE)
  #print(motif2)
  #If the k-mer is palindromic, consider the forward orientation even if the reverse is better
  
  
  if(class(motif1)!="try-error"){
    experiments_table$motif1[i]=motif1$motif
    experiments_table$orientation1[i]=motif1$kmer_orientation
    experiments_table$start1[i]=motif1$kmer_start
    experiments_table$end1[i]=motif1$kmer_end
  }
  
 
 
  
  #motif1=c(motif1,kmer_position_and_orientation(motif1, kmer1_freq_table, monomers))
  #motif2=c(motif2,kmer_position_and_orientation(motif2, kmer2_freq_table, monomers))
  
  if(class(motif2)!="try-error"){
    experiments_table$motif2[i]=motif2$motif
    experiments_table$orientation2[i]=motif2$kmer_orientation
    experiments_table$start2[i]=motif2$kmer_start
    experiments_table$end2[i]=motif2$kmer_end
  }
  
  i=i+1

}


write.table(experiments_table, file="~/projects/TFBS/PWMs_final_version2.2/kmer_frequencies_experiments_table_better.tsv", 
sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

table(experiments_table$orientation1)
#Forward Reverse 
#3    1312      79 

table(experiments_table$orientation2)
#Forward Reverse 
#1315      79


#experiments_table=read.table("~/projects/TFBS/PWMs_final_version2.2/kmer_frequencies_experiments_table_better.tsv", sep="\t", header=TRUE)

#which TFs fail

test=experiments_table %>% filter(motif1=="" | motif2=="" )

test=experiments_table %>% filter(motif1=="")

test=experiments_table %>% filter(motif2=="" )



#TF1s without motif
TF1s_without_motif=experiments_table %>% filter(motif1=="") %>% select(TF1) %>% unique() %>% pull(TF1)

TF1s_without_motif[!TF1s_without_motif %in% monomers$symbol]
#  "ZBTB33" "PBX2"   "gntR" 

TF1s_without_motif[!tolower(TF1s_without_motif) %in% tolower(monomers$symbol)]

#TF1s without motif
TF2s_without_motif=experiments_table %>% filter(motif2=="") %>% select(TF2) %>% unique() %>% pull(TF2)

TF2s_without_motif[!TF2s_without_motif %in% monomers$symbol]

TF2s_without_motif[!tolower(TF2s_without_motif) %in% monomers$symbol]

TF2s_without_motif[!tolower(TF2s_without_motif) %in% tolower(monomers$symbol)]

#These motifs need monomer
#"BHLHB8" CATATG BHLHA15
#"PROX2"  GGCGTC PROX1
#"NKX2-4" CACTTA NKX2-3
#"DBX2"   TAATTA DLX2
#"BHLHB5" CATATG BHLHA15
#"ZBTB33" ACATAA ZBTB33 ?
#"PBX2"   AATCAA PBX2 ?
#"FOXI2" TAAACA FOXI1
#"gntR" GTAACA gntR ?
#"NKX2-6" CACTTA NKX2-3

alternatives=list()
alternatives[["BHLHB8"]]="BHLHA15"
alternatives[["PROX2"]]="PROX1"
alternatives[["NKX2-4"]]="NKX2-3"
alternatives[["DBX2"]]=" DLX2"
alternatives[["BHLHB5"]]="BHLHA15"
alternatives[["ZBTB33"]]="ZBTB33" #?
alternatives[["PBX2"]]="PBX2" #?
alternatives[["FOXI2"]]="FOXI1"
alternatives[["gntR"]]="gntR" #?
alternatives[["NKX2-6"]]="NKX2-3"

missing_monomer=c("BHLHB8", "PROX2", "NKX2-4", "DBX2", "BHLHB5", "ZBTB33", "PBX2", "FOXI2", "gntR", "NKX2-6")

experiments_table %>% filter(TF1 %in% missing_monomer[1]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[2]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[3]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[4]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[5]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[6]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[7]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[8]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF1 %in% missing_monomer[9]) %>% select(kmer1) %>%unique()
experiments_table %>% filter(TF2 %in% missing_monomer[10]) %>% select(kmer2) %>%unique()


#list the 6-mers of each TF, are they unique for all

TFs=unique(c(experiments_table %>% pull(TF1), experiments_table %>% pull(TF2)))

TF_kmers=data.frame(TF=character(), kmer=character())

for(TF in TFs){
  #TF=TFs[1]
  #These seem to be unique
  #print(length(unique(c(experiments_table %>% filter(TF1==TF) %>% pull(kmer1), experiments_table %>% filter(TF2==TF) %>% pull(kmer2)))))
  kmer=unique(c(experiments_table %>% filter(TF1==TF) %>% pull(kmer1), experiments_table %>% filter(TF2==TF) %>% pull(kmer2)) )
  TF_kmers=rbind(TF_kmers, data.frame(TF=TF, kmer=kmer))

}

write.table(TF_kmers, file="~/projects/TFBS/PWMs_final_version2.2/TF_kmers_from_kmer_frequency_table.tsv", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
