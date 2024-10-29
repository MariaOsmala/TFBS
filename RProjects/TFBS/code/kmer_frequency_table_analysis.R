library("tidyverse")
library(readr)

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

length(unique(core_motifs$HGNC)) #404

TFs=names(which(table(core_motifs$HGNC)>1 ))
for(TF in TFs){
  #print(TF)
  #print(core_motifs %>% filter(HGNC==TF) %>% select(TF1.core.binding.site.1) %>% unique())
  if(nrow(core_motifs %>% filter(HGNC==TF) %>% select(TF1.core.binding.site.1) %>% unique())!=1){
    print(TF)
  }
}

#The 6-mer is different for the same TF depending on whether the TF is bait or prey

Prey_core_motifs %>% filter(HGNC=="HOXD12")
Bait_core_motifs %>% filter(HGNC=="HOXD12")


Prey_core_motifs %>% filter(HGNC=="SNAI1")
Bait_core_motifs %>% filter(HGNC=="SNAI1")

#Seems that in this table, if there are 4 orientations, the order is
#HT
#TH
#TT
#HH

#If there are two orientations, they are
#TH
#HT

#If there is only one orientation, it is
#HT

frequency_table <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/6mer_counts_fromDrive_28082024.csv", 
                                                                          delim = ";", escape_double = FALSE, trim_ws = TRUE) #3914 rows

frequency_table=frequency_table %>% unite(ID, c(TF1_TF2, Barcode, Batch),sep="_", remove=FALSE) #  "HOXC11_ERG"     "HOXC11_MEIS3"   "HOXC11_TEAD4"   "HOXC11_TGIF2LX" "NR4A2_FOXD2"    "TEAD4_CLOCK"

length(unique(frequency_table$ID)) #1394

names(table(frequency_table$ID)[which(table(frequency_table$ID)>4)]) #None

experiments=unique(frequency_table$ID)

#Heatmap for ERG_TBX6

count_matrix=frequency_table %>% filter(TF1_TF2=="ERG_TBX6") 
count_matrix=as.matrix(count_matrix[,-c(1:5)])
rownames(count_matrix)=frequency_table %>% filter(TF1_TF2=="ERG_TBX6") %>% select(kmer1.kmer2) %>% pull()
count_matrix=count_matrix/max(count_matrix)

colnames(count_matrix)=gsub("Count Gap ", "", colnames(count_matrix))

library(ComplexHeatmap)
library(circlize)

Heatmap(
  count_matrix[,1:15],
  name = "6-mer counts divided by max",            # Legend title
  col = colorRamp2(c( 0, 1), c("white", "red")),               # Custom color palette
  cluster_rows = FALSE,        # Do not cluster rows
  cluster_columns = FALSE,     # Do not cluster columns
  show_row_dend = FALSE,       # Hide row dendrogram
  show_column_dend = FALSE,     # Hide column dendrogram
  column_title="Gap",
  row_title="6-mers",
  row_title_side="left",
  column_title_side="top",
  row_names_side="left",
  column_names_side="top",
  column_names_rot=0
)


source("code/half_site_recognition_functions.R", echo=FALSE)

#Need halfsite for OTP, ERG

metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

#Metadata for artificial half sites, 7 motifs
metadata_halfsites <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

metadata_halfsites$new_representative="NO"
metadata_halfsites$filename=paste0("../../",metadata_halfsites$filename)
metadata_all=rbind(metadata[, colnames(metadata_halfsites)], metadata_halfsites)

#These were conserved in spacing motif analysis
#[1] "DLX1_MEIS3_NTGACANNNTAATTRN"            "DRGX_TEAD4_NWATTARCATWCCN"             
#[3] "ERG_TBX6_RSMGGAARYNNNRGGTGTNRN"         "GATA5_ONECUT1_NATCGATNNNNNNNNMGATAAN"  
#[5] "HOXB13_MEIS3_MATAAANNTGTCA"             "HOXC11_MEIS3_NYMRTAAANNTGTCAN"         
#[7] "HOXC11_MEIS3_NYMRTAAANNTGTCAN"          "LMX1A_BACH2_NRTGASTCANNNNNNNNNNTAATTRN"
#[9] "LMX1A_TGIF2LX_NYAATTANNTGTCAN"          "OLIG2_ISX_NTAATTRNNNNCAKMTGN"

#which(paste0(metadata$symbol, "_", metadata$seed)=="ERG_TBX6_RSMGGAARYNNNRGGTGTNRN") #2638


###################### Search TF pairs in the frequence table for which there is a spacing motif ##########################
#which kmer1.kmer2=="ACCGGA.TAAACA"

symbols1=strsplit(metadata %>% filter(experiment=="CAP-SELEX" & study=="fromYimeng" & type=="spacing") %>% select(symbol) %>% pull(symbol), "_")

padded_list <- lapply(symbols1, function(x) {
  length(x) <- 2
  return(x)
})

symbols1=do.call(rbind, padded_list)

paste0(symbols1[,1], "_",symbols1[,2])

head(apply(as.matrix(symbols1),1,paste0))

which(frequency_table %>% filter(kmer1.kmer2=="ACCGGA.TAAACA") %>%select(TF1_TF2) %>% pull(TF1_TF2) %in% paste0(symbols1[,1], "_",symbols1[,2]) )
which(frequency_table %>% filter(kmer1.kmer2=="ACCGGA.TAAACA") %>%select(TF1_TF2) %>% pull(TF1_TF2) %in% paste0(symbols1[,2], "_",symbols1[,1]) )
frequency_table %>% filter(kmer1.kmer2=="ACCGGA.TAAACA") %>%select(TF1_TF2) %>% pull(TF1_TF2)

test=frequency_table[which(frequency_table %>%select(TF1_TF2) %>% pull(TF1_TF2) %in% paste0(symbols1[,1], "_",symbols1[,2]) ),]

test2=frequency_table[which(frequency_table %>%select(TF1_TF2) %>% pull(TF1_TF2) %in% paste0(symbols1[,2], "_",symbols1[,1]) ),]

######################################################################################


monomers=metadata_all %>% filter(experiment=="HT-SELEX") # 1745

#problem is k-mers in different order or some other error
experiments_table=data.frame(experiment=experiments, number_of_orientations=NA, in_heterodimers=NA, problem="", kmer1="", motif1="", orientation1="", position1="", 
                             kmer2="", motif2="", orientation2="", position2="")






i=1
for(ex in experiments){
  #ex=experiments[1]
  print(ex)
  data=frequency_table %>% filter(ID==ex) %>% select(-ID)
  
  print(paste0( nrow(data), "_orientations"))
  
  experiments_table[i, "number_of_orientations"]=as.numeric(nrow(data))
  
  TF1=unlist(strsplit(data$TF1_TF2[1], "_"))[1]
  TF2=unlist(strsplit(data$TF1_TF2[1], "_"))[2]
  
  kmer1_freq_table=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][1] #"TAAACA"
  kmer2_freq_table=strsplit(data$kmer1.kmer2[1], "[.]")[[1]][2] #"TGCGGG"
  
  #In heterodimers?
  experiments_table$in_heterodimers[i]=data$TF1_TF2[1] %in% Prey_core_motifs_heterodimers$HGNC
  
  #In baits and prey
  TF1 %in% Prey_core_motifs$HGNC #TRUE
  kmer1_core_Prey=as.character(Prey_core_motifs %>% filter(HGNC==TF1) %>% select(TF1.core.binding.site.1)) #"TAAACA"
  TF1 %in% Bait_core_motifs$HGNC #FALSE
  kmer1_core_Bait=as.character(Bait_core_motifs %>% filter(HGNC==TF1) %>% select(TF1.core.binding.site.1))
  
  TF2 %in% Prey_core_motifs$HGNC #TRUE
  kmer2_core_Prey=as.character(Prey_core_motifs %>% filter(HGNC==TF2) %>% select(TF1.core.binding.site.1)) #"TGCGGG"
  TF2 %in% Bait_core_motifs$HGNC #TRUE
  kmer2_core_Bait=as.character(Bait_core_motifs %>% filter(HGNC==TF2) %>% select(TF1.core.binding.site.1)) #"TGCGGG"
  
  
  
  
  
  
  #Check that the kmer is the same in the frequency table as in the core motifs
  
  if(kmer1_core_Prey=="character(0)" & kmer1_core_Bait!="character(0)"){
    kmer1_core_Prey=kmer1_core_Bait
  } else if(kmer1_core_Prey!="character(0)" & kmer1_core_Bait=="character(0)"){
    kmer1_core_Bait=kmer1_core_Prey
  }else if(kmer1_core_Prey=="character(0)" & kmer1_core_Bait=="character(0)"){
    print("The 6-mer of TF1 is not in the core motifs")
  }
  
  if(kmer2_core_Prey=="character(0)" & kmer2_core_Bait!="character(0)"){
    kmer2_core_Prey=kmer2_core_Bait
  } else if(kmer2_core_Prey!="character(0)" & kmer2_core_Bait=="character(0)"){
    kmer2_core_Bait=kmer2_core_Prey
  }else if(kmer2_core_Prey=="character(0)" & kmer2_core_Bait=="character(0)"){
    print("The 6-mer of TF2 is not in the core motifs")
  }
  
  if(kmer1_core_Prey==kmer1_freq_table && kmer1_core_Bait==kmer1_freq_table){
    print("kmers of TF1 match")
  }else{
    print("kmers of TF1 do not match")
  }
  
  if(kmer2_core_Prey==kmer2_freq_table && kmer2_core_Bait==kmer2_freq_table){
    print("kmers of TF2 match")
  }else{
    print("kmers of TF2 do not match")
  }
  
  #If k-mers do not match, are they in different order?
  
  if(!(kmer1_core_Prey==kmer1_freq_table && kmer1_core_Bait==kmer1_freq_table && kmer2_core_Prey==kmer2_freq_table && kmer2_core_Bait==kmer2_freq_table )){
    
    tmp=kmer1_freq_table
    kmer1_freq_table=kmer2_freq_table
    kmer2_freq_table=tmp
    if((kmer1_core_Prey==kmer1_freq_table && kmer1_core_Bait==kmer1_freq_table && kmer2_core_Prey==kmer2_freq_table && kmer2_core_Bait==kmer2_freq_table )){
      print("kmers in different order")
      experiments_table$problem[i]="kmers in different order"
    }else{
      print("some other error than k-mers in different order")
      tmp=kmer1_freq_table
      kmer1_freq_table=kmer2_freq_table
      kmer2_freq_table=tmp
      experiments_table$problem[i]="some other problem"
    }
  }
  
  if(kmer1_core_Prey==kmer1_freq_table && kmer1_core_Bait==kmer1_freq_table){
    print("kmers of TF1 match")
  }else{
    print("kmers of TF1 still do not match, use the kmer1_freq_table")
    
  }
  
  if(kmer2_core_Prey==kmer2_freq_table && kmer2_core_Bait==kmer2_freq_table){
    print("kmers of TF2 match")
  }else{
    print("kmers of TF2 still do not match, use the kmer1_freq_table")
  }
  
  
  #Head-to-tail "TAAACA.TGCGGG" F1 F2
  #Tail-to-head "TGCGGG.TAAACA" F2 F1
  #Tail-to-tail "TGTTTA.TGCGGG" RC1 F2
  #Head-to-head "TAAACA.CCCGCA" F1 RC2
  
  #Find the motifs. The motif should match with the k-mer
  
  experiments_table$kmer1[i]=kmer1_freq_table
  experiments_table$kmer2[i]=kmer2_freq_table
  
  motif1=try(find_monomer_match_kmer(monomer=TF1, monomers, kmer1_freq_table), silent=TRUE)
  #print(motif1)
  motif2=try(find_monomer_match_kmer(monomer=TF2, monomers, kmer2_freq_table), silent=TRUE)
  #print(motif2)
  #If the k-mer is palindromic, consider the forward orientation even if the reverse is better
  
  if( class(motif1)=="try-error" || class(motif2)=="try-error"){
    i=i+1
    next
    
  }
  

  experiments_table$motif1[i]=motif1$motif
  experiments_table$motif2[i]=motif2$motif
 
  
  motif1=c(motif1,kmer_position_and_orientation(motif1, kmer1_freq_table, monomers))
  motif2=c(motif2,kmer_position_and_orientation(motif2, kmer2_freq_table, monomers))
  
  
  experiments_table$orientation1[i]=motif1$orientation
  experiments_table$position1[i]=motif1$position
  
  
  
  experiments_table$orientation2[i]=motif2$orientation
  experiments_table$position2[i]=motif2$position
  
  
  
  if(motif1$orientation=="Reverse"){
    print("kmer1 reverse orientation")
  }
  
  if(motif2$orientation=="Reverse"){
    print("kmer2 reverse orientation")
  }
  
  i=i+1

}


write.table(experiments_table, file="~/projects/TFBS/PWMs_final_version2.2/kmer_frequencies_experiments_table.tsv", 
sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


experiments_table=read.table("~/projects/TFBS/PWMs_final_version2.2/kmer_frequencies_experiments_table.tsv", sep="\t", header=TRUE)

experiments_table$number_of_orientations %>% table()
#1   2   4 
#94 690 610 


experiments_table$problem %>% table()
#no problem                                                   kmers in different order       some other problem 
#544 (these are all with 4 orientations)                      768                            82 

experiments_table %>% filter(problem=="") %>% select(number_of_orientations) %>% table() 
#4 
#544 

experiments_table %>% filter(problem=="kmers in different order") %>% select(number_of_orientations) %>% table() 
#1   2 
#94 674 
  
experiments_table %>% filter(problem=="some other problem") %>% select(number_of_orientations) %>% table() 
#2  4 
#16 66 