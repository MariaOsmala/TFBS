library("tidyverse")
library(readr)
library("readxl")
source("code/half_site_recognition_functions.R", echo=FALSE)

#Which package has read_excel function 

#install.packages("readxl")

Jolma2015_6mers <- read_excel("~/projects/TFBS/PWMs_final_version2.2/6mers_from_Jolma2015/6mers.xlsx")
#Jolma2015_6mers$`6-mer`


Bait_core_motifs=read.table("../../PWMs_final_version2.2/fromYimeng/core_6mers/Bait.core_motif", sep="\t", header=TRUE) #160

length(unique(Bait_core_motifs$HGNC)) #158
which(table(Bait_core_motifs$HGNC)>1 )

#GMEB2 TP73L These are twice in the table but with the same 6-mer
#51   147 

Bait_core_motifs %>% filter(HGNC=="GMEB2") 
Bait_core_motifs %>% filter(HGNC=="TP73L")

#Remove duplicate rows from the table
Bait_core_motifs=Bait_core_motifs %>% distinct() #158

Prey_core_motifs=read.table("../../PWMs_final_version2.2/fromYimeng/core_6mers/Prey.core_motif", sep="\t", header=TRUE) #368

length(unique(Prey_core_motifs$HGNC)) #364
which(table(Prey_core_motifs$HGNC)>1 )

Prey_core_motifs %>% filter(HGNC=="BCL11B") 
Prey_core_motifs %>% filter(HGNC=="GFI1B")
Prey_core_motifs %>% filter(HGNC=="HOXC11")
Prey_core_motifs %>% filter(HGNC=="NR4A2")

#Remove duplicate rows from the table
Prey_core_motifs=Prey_core_motifs %>% distinct() #364

#Separate heterodimers 

Prey_core_motifs_heterodimers=Prey_core_motifs[grep("_",Prey_core_motifs$HGNC),]
Prey_core_motifs=Prey_core_motifs[-grep("_",Prey_core_motifs$HGNC),]
#remove column 
Prey_core_motifs=Prey_core_motifs %>% select( -TF2.core.binding.site)

Prey_core_motifs_heterodimers=cbind(Prey_core_motifs_heterodimers, do.call(rbind, strsplit(Prey_core_motifs_heterodimers$HGNC, "_")))

names(Bait_core_motifs)=c("HGNC", "TF1.core.binding.site.1")
Bait_core_motifs$TF1.core.binding.site.2=NA

core_motifs=rbind(Bait_core_motifs, Prey_core_motifs)


#Generated in kmer_frequency_table_analysis_better.R
TF_kmers=read.table(file="~/projects/TFBS/PWMs_final_version2.2/TF_kmers_from_kmer_frequency_table.tsv",header=TRUE, sep="\t")


metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

capselex=metadata %>% filter(experiment=="CAP-SELEX") #1898


symbols=as.data.frame(do.call(rbind, strsplit(capselex$symbol, "_")))
names(symbols)=c("TF1", "TF2")

capselex <- capselex %>%
  mutate(TF1 = symbols$TF1) %>%  
  relocate(TF1, .after = symbol)     

capselex <- capselex %>%
  mutate(TF2 = symbols$TF2) %>%  
  relocate(TF2, .after = TF1)    

families=as.data.frame(do.call(rbind, strsplit(capselex$Lambert2018_families, "_")))
names(families)=c("TF1", "TF2")
capselex <- capselex %>%
  mutate(TF1_family =families$TF1) %>%  
  relocate(TF1_family, .after = TF2)     

capselex <- capselex %>%
  mutate(TF2_family = families$TF2) %>%  
  relocate(TF2_family, .after = TF1_family)    


capselex<- capselex %>%
  mutate(kmer1="") %>%
  relocate(kmer1, .after = TF1_family)  

capselex<- capselex %>%
  mutate(kmer2="") %>%
  relocate(kmer2, .after = kmer1)  

capselex<- capselex %>%
  mutate(order_orientation="") %>%
  relocate(order_orientation, .after = kmer2)  


capselex<- capselex %>%
  mutate(kmer_spacing="") %>%
  relocate(kmer_spacing, .after = order_orientation)

#Check that these go right
# KLF13_FLI1
# KLF14_FLI1

#which(capselex$symbol=="KLF13_FLI1") #1226 1227
#which(capselex$symbol=="KLF14_FLI1") #1228 1229
#which(capselex$ID=="KLF13_FLI1_TGATGA40NCATT_YYIII_NRNCACGCCCCTTCCGNYN_m2_c3b0")
#which(capselex$ID=="POU3F2_MEIS3_TTGTCC40NTAG_YRIIII_NTGACANTAATKWN_m1_c3b0")
#i=1879

missing_kmer1=c()
missing_kmer2=c()

for(i in 1:nrow(capselex)){
  #i=1
  print(i)
  
  TF1=capselex$TF1[i]
  TF2=capselex$TF2[i]
  
  kmer1=TF_kmers %>% filter(TF==TF1) %>% pull(kmer)
  kmer2=TF_kmers %>% filter(TF==TF2) %>% pull(kmer)
  
  if(length(kmer1)==0){
    
    if(capselex$study[i]=="fromYimeng"){
      print("no kmer1 in fromYimeng")
      #These are all composite motifs for which the probably is no 6-mer table
      missing_kmer1=c(missing_kmer1, i)
      kmer1=unique(core_motifs %>% filter(HGNC==TF1) %>% pull(TF1.core.binding.site.1)) 
    }
    if(capselex$study[i]=="Jolma2015"){
        kmer1=unique(Jolma2015_6mers %>% filter(`HNGC-name`==TF1) %>% pull(`6-mer`))
    }
    
  }
  
  if(length(kmer2)==0){
    
    if(capselex$study[i]=="fromYimeng"){
      print("no kmer2 in fromYimeng")
      missing_kmer2=c(missing_kmer2, i)
      kmer2=unique(core_motifs %>% filter(HGNC==TF2) %>% pull(TF1.core.binding.site.1)) 
  
    }
    if(capselex$study[i]=="Jolma2015"){
      kmer2=unique(Jolma2015_6mers %>% filter(`HNGC-name`==TF2) %>% pull(`6-mer`))
    }
    

    
  }
  
  capselex$kmer1[i]=kmer1
  capselex$kmer2[i]=kmer2
  
  kmer1_positions=kmer_position_in_motif(motif_filename=capselex$filename[i], kmer=kmer1)
  kmer2_positions=kmer_position_in_motif(motif_filename=capselex$filename[i], kmer2)
  
  order_orientation=""
  spacing=NA
  
  if(kmer1_positions$kmer_start<kmer2_positions$kmer_start){
    #TF1_TF2
    
    if(kmer1_positions$kmer_orientation=="Forward" & kmer2_positions$kmer_orientation=="Forward"){
      order_orientation="HT"
      
    }else if(kmer1_positions$kmer_orientation=="Reverse" & kmer2_positions$kmer_orientation=="Reverse"){
      order_orientation="TH"
    }else if(kmer1_positions$kmer_orientation=="Forward" & kmer2_positions$kmer_orientation=="Reverse"){
      order_orientation="HH"
    }else{
      order_orientation="TT"
      #kmer1_positions$kmer_orientation=="Reverse" & kmer2_positions$kmer_orientation=="Forward"
    }
    
    spacing=kmer2_positions$kmer_start-kmer1_positions$kmer_end-1
    
    
  }else{
    #TF2_TF1
    
    if(kmer1_positions$kmer_orientation=="Reverse" & kmer2_positions$kmer_orientation=="Reverse"){
      order_orientation="HT"
      
    }else if(kmer1_positions$kmer_orientation=="Forward" & kmer2_positions$kmer_orientation=="Forward"){
      order_orientation="TH"
    }else if(kmer1_positions$kmer_orientation=="Reverse" & kmer2_positions$kmer_orientation=="Forward"){
      order_orientation="HH"
    }else{
      order_orientation="TT"
      #kmer1_positions$kmer_orientation=="Forward" & kmer2_positions$kmer_orientation=="Reverse"
    }
    
    
    spacing=kmer1_positions$kmer_start-kmer2_positions$kmer_end-1
    
    
  }
  
  
  capselex$order_orientation[i]=order_orientation
  capselex$kmer_spacing[i]=spacing
  
  
    
}
  

write.table(capselex, file="~/projects/TFBS/PWMs_final_version2.2/capselex_kmer_order_orientation_spacing.tsv", 
sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


ETV2_FOXI1=capselex %>%filter(symbol %in% c("FOXI1_ETV2", "ETV2_FOXI1")) %>%
  mutate(full_type= str_c(order_orientation, kmer_spacing, sep = ""))  %>% select(c(ID,full_type, type, new_representative))

#What represents FOXI1_ETV2 motifs?

represented_by_representatives <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/represented_by_representatives_version2.2.RDS")

ETV2_FOXI1$representative_motif=""
ETV2_FOXI1$representative_motif[1]=names(represented_by_representatives)[sapply(represented_by_representatives, function(x) ETV2_FOXI1$ID[1] %in% x)]
ETV2_FOXI1$representative_motif[2]=names(represented_by_representatives)[sapply(represented_by_representatives, function(x) ETV2_FOXI1$ID[2] %in% x)]
ETV2_FOXI1$representative_motif[3]=names(represented_by_representatives)[sapply(represented_by_representatives, function(x) ETV2_FOXI1$ID[3] %in% x)]
ETV2_FOXI1$representative_motif[4]=names(represented_by_representatives)[sapply(represented_by_representatives, function(x) ETV2_FOXI1$ID[4] %in% x)]
ETV2_FOXI1$representative_motif[5]=names(represented_by_representatives)[sapply(represented_by_representatives, function(x) ETV2_FOXI1$ID[5] %in% x)]

#Onko näille kmer-frequency tablea?

kmer_frequencies_experiments_table=read.table("~/projects/TFBS/PWMs_final_version2.2/kmer_frequencies_experiments_table_better.tsv", 
                                              sep="\t", header=TRUE)

which((capselex %>% filter(ID %in% ETV2_FOXI1$representative_motif) %>% pull(symbol)) %in% (kmer_frequencies_experiments_table %>% mutate(symbol= str_c(TF1, TF2, sep = "_")) %>%pull(symbol)))
which((capselex %>% filter(ID %in% ETV2_FOXI1$representative_motif) %>% pull(symbol)) %in% (kmer_frequencies_experiments_table %>% mutate(symbol= str_c(TF2, TF1, sep = "_")) %>%pull(symbol)) )

(capselex %>% filter(ID %in% ETV2_FOXI1$representative_motif) %>% pull(ID))[1] #representative

#is the kmer-frequncy table present for the motifs represented by representatives

nonreps=unlist(represented_by_representatives[ETV2_FOXI1$representative_motif])
names(nonreps)=NULL

nonreps=unique(nonreps)
which((capselex %>% filter(ID %in% nonreps) %>% pull(symbol)) %in% (kmer_frequencies_experiments_table %>% mutate(symbol= str_c(TF1, TF2, sep = "_")) %>%pull(symbol)))
which((capselex %>% filter(ID %in% nonreps) %>% pull(symbol)) %in% (kmer_frequencies_experiments_table %>% mutate(symbol= str_c(TF2, TF1, sep = "_")) %>%pull(symbol)) )

(capselex %>% filter(ID %in% nonreps) %>% pull(ID))[c(1,2)] #non-representative


"ELK1_FOXI1_CAP-SELEX_TAGTCA40NTCT_AAA_RSCGGATGTKKN_1_3"  
"ELK1_FOXI1_CAP-SELEX_TAGTCA40NTCT_AAA_WGTTKMCGGAWRTN_1_3"
"ELK1_FOXI1_CAP-SELEX_TAGTCA40NTCT_AAA_RSCGGAANNRTMAAYAN_1_3"

tmp=capselex %>% filter(symbol=="ELK1_FOXI1")

#which( (kmer_frequencies_experiments_table %>% mutate(symbol= str_c(TF1, TF2, sep = "_")) %>%pull(symbol) )=="FOXI1_ELK1" ) #664
