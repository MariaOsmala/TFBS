

rm(list=ls())


find_common_prefix <- function(strings) {
  # If the vector is empty, return an empty string
  if (length(strings) == 0) return("")
  
  # Find the minimum length of strings
  min_length <- min(nchar(strings))
  
  # Initialize common prefix as an empty string
  common_prefix <- ""
  
  # Iterate through each character position
  for (i in 1:min_length) {
    # Get the i-th character from each string
    chars <- substr(strings, i, i)
    
    # If not all characters are the same, break
    if (length(unique(chars)) > 1) {
      break
    }
    
    # Append the common character to the common prefix
    common_prefix <- paste0(common_prefix, chars[1])
  }
  
  return(common_prefix)
}

library(readr)
library(tidyverse)

#TFBS_path='/Users/osmalama/projects/TFBS/'
# These were able to be computed only part of the motif pairs ~3 million, there should be 16 million (both directions)
tomtom = read_tsv('/Users/osmalama/projects/TFBS/Results/tomtom_final/tomtom.all.txt') #2_722_099 rows
#remove three last rows
tomtom <- tomtom %>% slice(1:(n()-3)) #2044862


metadata <- read_delim("~/projects/TFBS/PWMs_final/metadata_representatives.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

#monomers

monomers=metadata[metadata$experiment!="CAP-SELEX",]

monomers_representative=monomers[which(monomers$new_representative=="YES"),]

#How to check which motif is heterodimeric and which is composite?

representatives=metadata[which(metadata$new_representative=="YES"),]

#represented_by_representatives <- read_table("~/projects/TFBS/PWMs_final/sstat_represented_by_representatives.tsv")

represented_by_representatives <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/represented_by_representatives.RDS")

#which are heterodimeric or composite motifs

capselex=representatives[representatives$experiment=="CAP-SELEX",]

capselex$type[grep("composite",capselex$ID)]="composite"
capselex$type[grep("spacing",capselex$ID)]="spacing"
capselex$type[grep("20230420",capselex$type)]="spacing"

capselex$first.monomer=""
capselex$second.monomer=""

capselex$First_representative_bool=NA
capselex$First_tomtom_reverse_order=NA
capselex$First_Optimal_offset=NA
capselex$First_Optimal_Overlap=NA
capselex$First_Optimal_Orientation=""

capselex$Second_representative_bool=NA
capselex$Second_tomtom_reverse_order=NA
capselex$Second_Optimal_offset=NA
capselex$Second_Optimal_Overlap=NA
capselex$Second_Optimal_Orientation=""

find_representative_monomer <- function(monomer, monomers, represented_by_representatives, metadata) {
  representative_bool=FALSE
  mono1=monomers[grep(monomer,monomers$ID),]
  
  #There can be representative dimers, remove those
  if(length(which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                     "putatively multimeric")))!=0){
    mono1=mono1[-which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                          "putatively multimeric")),]
  }
  
  
  #Exclude METHYL-SELEX
  if(length(which(mono1$experiment=="Methyl-HT-SELEX"))!=0){
    mono1=mono1[-which(mono1$experiment=="Methyl-HT-SELEX"),]
  }
  
  if( length(which(mono1$new_representative=="YES"))!=0 ){
    mono1=mono1[which(mono1$new_representative=="YES"),]
    representative_bool=TRUE
  }
  common_prefix=find_common_prefix(mono1$ID[grep(monomer,mono1$ID)])
  if(representative_bool==FALSE){
    indices <- which(sapply(represented_by_representatives, function(x) any(grepl(common_prefix, x))))
    #can be many, Let's take the shortes, or the one with highest information content?
    mono_representative=names(represented_by_representatives)[indices]
    motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) #Atoh1
    if(length(motif1)>1){
      motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>%filter(IC==max(IC))%>% pull(ID)
    }
    
  }else{
    motif1=mono1 %>%filter(length==min(length)) %>% pull(ID) 
    if(length(motif1)>1){
      motif1=mono1 %>%filter(length==min(length)) %>%filter(IC==max(IC)) %>% pull(ID)
    }
  }
  
  list(motif=motif1, representative_bool=representative_bool)
}


for(i in 1:nrow(capselex)){
  print(i)
  dimer_ID=capselex$ID[i]
  #monomers
  monomers_symbols=unlist(strsplit(capselex$symbol[i],"_"))
  
  #First monomer, is any of these representative
  
  motif1=find_representative_monomer(monomer=monomers_symbols[1], monomers, represented_by_representatives, metadata)
  motif2=find_representative_monomer(monomer=monomers_symbols[2], monomers, represented_by_representatives, metadata)
  
  tomtom_reverse_order=FALSE
  indices <- which(
    (tomtom$Query_ID==motif1$motif) & (tomtom$Target_ID==dimer_ID) & 
      tomtom$Query_ID != tomtom$Target_ID
  )
  
  if(length(indices)==0){
    tomtom_reverse_order=TRUE
    indices <- which(
      (tomtom$Target_ID==motif1$motif) & (tomtom$Query_ID==dimer_ID) & 
        tomtom$Query_ID != tomtom$Target_ID
    )
  }
  
  capselex$first.monomer[i]=motif1$motif
  capselex$second.monomer[i]=motif2$motif
  
  capselex$First_representative_bool[i]=motif1$representative_bool
  capselex$First_tomtom_reverse_order[i]=tomtom_reverse_order
  if(length(indices)!=0){
    capselex[i, c("First_Optimal_offset","First_Optimal_Overlap","First_Optimal_Orientation")]=tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")]
  }
  tomtom_reverse_order=FALSE
  indices <- which(
    (tomtom$Query_ID==motif2$motif) & (tomtom$Target_ID==dimer_ID) & 
      tomtom$Query_ID != tomtom$Target_ID
  )
  
  if(length(indices)==0){
    tomtom_reverse_order=TRUE
    indices <- which(
      (tomtom$Target_ID==motif2$motif) & (tomtom$Query_ID==dimer_ID) & 
        tomtom$Query_ID != tomtom$Target_ID
    )
  }
  
  capselex$Second_representative_bool[i]=motif2$representative_bool
  capselex$Second_tomtom_reverse_order[i]=tomtom_reverse_order
  if(length(indices)!=0){
    capselex[i, c("Second_Optimal_offset","Second_Optimal_Overlap","Second_Optimal_Orientation")]=tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")]
  }
  
}

#match for both motifs

monomer_start_end <- function(tomtom_reverse_order, optimal_offset, optimal_orientation, optimal_overlap, dimer_length) {
  
  #tomtom_reverse_order=capselex$First_tomtom_reverse_order[i]
  #optimal_offset=capselex$First_Optimal_offset[i]
  #optimal_orientation=capselex$First_Optimal_Orientation[i]
  #optimal_overlap=capselex$First_Optimal_Overlap
  #dimer_length=capselex$length[i]
  
  #tomtom_reverse_order=capselex$Second_tomtom_reverse_order[i]
  #optimal_offset=capselex$Second_Optimal_offset[i]
  #optimal_orientation=capselex$Second_Optimal_Orientation[i]
  #optimal_overlap=capselex$Second_Optimal_Overlap[i]
  #dimer_length=capselex$length[i]
  
  
  if(tomtom_reverse_order==FALSE){
    offset=optimal_offset
  }else{
    offset=-optimal_offset
  }
  if(optimal_orientation=="-"){
    if(offset<0){
      offset=0
    }
    start=dimer_length-offset-optimal_overlap+1
    end=dimer_length-offset
  }else{
    if(offset<0){
      offset=0
    }
    start=1+offset
    end=offset+optimal_overlap
  }
  #The start needs to be min 1 and end needs to be max the dimer_length of the dimer
  
  start=max(1, start)
  end=min(dimer_length, end)
  
  list(start=start, end=end)
}


capselex$type_computation=""
capselex$overlap_length=NA
capselex$overlap_start=NA
capselex$overlap_end=NA

capselex$gap_length=NA
capselex$gap_start=NA
capselex$gap_end=NA



for (i in which( paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation) %in% c("- -", "- +", "+ +", "+ -")) ){
  print(i)
  #i=which( paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation) %in% c("- -", "- +", "+ +", "+ -"))[1]
  #"first.monomer"               "second.monomer"              "First_representative_bool"   "First_tomtom_reverse_order"  "First_Optimal_offset"       
  #"First_Optimal_Overlap"       "First_Optimal_Orientation"   "Second_representative_bool"  "Second_tomtom_reverse_order" "Second_Optimal_offset"      
  #"Second_Optimal_Overlap"      "Second_Optimal_Orientation" 
  
  first_se=monomer_start_end( tomtom_reverse_order=capselex$First_tomtom_reverse_order[i], 
                     optimal_offset=capselex$First_Optimal_offset[i], 
                     optimal_orientation=capselex$First_Optimal_Orientation[i], 
                     optimal_overlap=capselex$First_Optimal_Overlap[i], 
                     dimer_length=capselex$length[i])
  
  
  second_se=monomer_start_end( tomtom_reverse_order=capselex$Second_tomtom_reverse_order[i], 
                              optimal_offset=capselex$Second_Optimal_offset[i], 
                              optimal_orientation=capselex$Second_Optimal_Orientation[i], 
                              optimal_overlap=capselex$Second_Optimal_Overlap[i], 
                              dimer_length=capselex$length[i])
  
  
  #cap length
  #Left and rightmost
  if(first_se$start< second_se$start){
    
    gap=second_se$start-first_se$end #0
    
    if(gap>0){
      capselex$type_computation[i]="spacing"
      capselex$gap_length[i]=second_se$start-first_se$end-1
      capselex$gap_start[i]=first_se$end+1
      capselex$gap_end[i]=second_se$start-1
    }else{
      capselex$type_computation[i]="composite"
      #The size of the overlapping part is abs(-3)+1=3
      capselex$overlap_length[i]=abs(gap)+1
      capselex$overlap_start[i]=second_se$start
      capselex$overlap_end[i]=first_se$end
      
    }
      
   }else{
     print("reverse")
     gap=first_se$start-second_se$end #0
     
     if(gap>0){
       capselex$type_computation[i]="spacing"
       capselex$gap_length[i]=first_se$start-second_se$end-1
       capselex$gap_start[i]=second_se$end+1
       capselex$gap_end[i]=first_se$start-1
     }else{
       capselex$type_computation[i]="composite"
       #The size of the overlapping part is abs(-3)+1=3
       capselex$overlap_length[i]=abs(gap)+1
       capselex$overlap_start[i]=first_se$start
       capselex$overlap_end[i]=second_se$end
       
     }
     
     
    
  }
  
    
  
}

#E2F3_HT-SELEX_TTACAT20NACA_AF_NNAAATGGCGCCAAAANN_2_4  -6                    12                         +  FALSE
#ELK1_HT-SELEX_TCGGAA20NAGT_AG_ACCGGAAGTN_1_2       8                      9  +                            FALSE
#In this case the coordinates are: 
#motif1 (E2F3): start: 1+ -6 =-5 (->0) end: 0(start)+12(overlap)=12
#motif2 (ELK1): start: 1+8=9 end: 8(offset)+9 (overlap)=17
#cap size: start(E2F3)-end(ELK1)= 11 - 8-1=2  
#cap size: start(ELK1)-end(E2F3)= 9-12=-3
#The size of the overlapping part is abs(-3)+1=3
#Coordinates of the overlapping part: start(ELK1)-end(E2F3)

#ALX4_HT-SELEX_TCTATT40NCAT_KW_CYAATTAN_1_3 12 8 - FALSE
#TBX3_HT-SELEX_TAAGCC40NAGT_AR_AGGTGTNR_1_4 0 8 + FALSE
#ALX4 start: 25-12-8+1=6    length(dimer)25-12=13 (orientation is negative)
#TBX3 start: 0+1 - 8

#TEAD3_HT-SELEX_TATAAG40NAAG_AI_RMATWCCA_1_3  0 8 - FALSE
#CEBPG_HT-SELEX_TAAAAT20NCG_AC_NTTRCGCAAY_1_2 12 10 -                    FALSE                      FALSE                    0
#TEAD3: start=length(dimer)22-0-8+1=15 end=length(dimer)-0=22
#CEBPG: start=length(dimer)22-12-10+1=0 end=length(dimer)-12=22-12=10




unique(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)) #"- -" "- +" "+ +" "+ -" 
table(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)) #"- -" "- +" "+ +" "+ -" 

length(which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="+ +")) #213
length(which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="- +")) #109
length(which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="+ -")) #114
length(which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="- -")) #32

plusplus=capselex[which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="+ +"),] #213
#Offsets can be positive or negative, both positive is clear

as.data.frame(capselex[capselex$ID=="HOXA10_TBX4_TGCGGT40NTCA_YWII_NAGGTGTNANRTCGTAAAN_m2_c3b0u_short_composite_new", 25:36])
#Hoxc10_HT-SELEX_TATAGA20NGAG_Z_GYAATWAAAN_1_4     10 9 + FALSE
#TBX3_HT-SELEX_TAAGCC40NAGT_AR_AGGTGTNR_1_4    1 8 + FALSE


#as.data.frame(capselex[capselex$ID=="ATOH1_ONECUT1_TCCATG40NCAC_YJII_CATATGNNNNNNATCGAT_m1_c3b0_short_spacing_new", 25:36])
#Atoh1_HT-SELEX_TAGTCC20NCTT_Z_RNCATATGNY_1_3 #10 8 - FALSE  
#ONECUT1_HT-SELEX_TAGCTC20NTCT_Y_NNAAAATCRATAWN_1_3     #7 11 + FALSE             



#First offset negative, other positive
as.data.frame(capselex[capselex$ID=="E2F1_ELK1_CAP-SELEX_TATGTA40NTGT_AX_SGCGCNNNNNCGGAAGN_1_2",25:36])
#E2F3_HT-SELEX_TTACAT20NACA_AF_NNAAATGGCGCCAAAANN_2_4  -6                    12                         +  FALSE
#ELK1_HT-SELEX_TCGGAA20NAGT_AG_ACCGGAAGTN_1_2       8                      9  +                            FALSE

#In this case the coordinates are: 
#motif1 (E2F3): start: 1+ -6 =-5 (->0) end: 0(start)+12(overlap)=12
#motif2 (ELK1): start: 1+8=9 end: 8(offset)+9 (overlap)=17
#cap position: end(ALX4)+1:start(Atoh1  )
#cap size: start(E2F3)-end(ELK1)= 11 - 8-1=2  
  
#cap size: start(ELK1)-end(E2F3)= 9-12=-3
#The size of the overlapping part is abs(-3)+1=3
#Coordinates of the overlapping part: start(ELK1)-end(E2F3)

#ALX4 start: 25-12-8+1=6    length(dimer)25-12=13 (orientation is negative)
#TBX3 start: 0+1 - 8
#TEAD3: start=length(dimer)22-0-8-1=15 end=length(dimer)-0=22
#CEBPG: start=length(dimer)22-12-10+1=0 end=length(dimer)-12=22-12=10
  
#First offset positive, second negative
as.data.frame(capselex[capselex$ID=="ELK1_FOXI1_CAP-SELEX_TAGTCA40NTCT_AAA_RSCGGAANRWMAACAN_1_3",25:36])
#ELK1_HT-SELEX_TCGGAA20NAGT_AG_ACCGGAAGTN_1_2 0 10 + #FALSE
#FOXA1_HT-SELEX_TTCTAA40NAAT_KN_TRNGTAAACA_1_3b1    -5-> 5 10 + TRUE, in this case offset needs to be of opposite sign
  

minusplus=capselex[which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="- +"),] 

as.data.frame(capselex[capselex$ID=="ALX4_TBX21_CAP-SELEX_TGGCAC40NAAG_AAD_RGGTGNTAATNNNNNNNNNCASYNN_1_3", 25:36])
#ALX4_HT-SELEX_TCTATT40NCAT_KW_CYAATTAN_1_3 12 8 - FALSE
#TBX3_HT-SELEX_TAAGCC40NAGT_AR_AGGTGTNR_1_4 0 8 + FALSE

#ALX4 start: 25-12-8+1=6    length(dimer)25-12=13 (orientation is negative)
#TBX3 start: 0+1 - 8
#TEAD3: start=length(dimer)22-0-8-1=15 end=length(dimer)-0=22
#CEBPG: start=length(dimer)22-12-10+1=0 end=length(dimer)-12=22-12=10

minusminus=capselex[which(paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation)=="- +"),] 
as.data.frame(capselex[capselex$ID=="TEAD4_CEBPB_CAP-SELEX_TGACTC40NTCA_AY_NTTRCGYAANNNNNNGGAATGY_1_3",25:36])
#TEAD3_HT-SELEX_TATAAG40NAAG_AI_RMATWCCA_1_3  0 8 - FALSE
#CEBPG_HT-SELEX_TAAAAT20NCG_AC_NTTRCGCAAY_1_2 12 10 -                    FALSE                      FALSE                    0

#TEAD3: start=length(dimer)22-0-8-1=15 end=length(dimer)-0=22
#CEBPG: start=length(dimer)22-12-10+1=0 end=length(dimer)-12=22-12=10

#cap size: start(TEAD3)-end(CEBPG)= 15-10=5

#TEAD3    start: 

#spacing motif, this is not representative
motif="OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0_short_20230420"

#monomeric motifs for OLIG2 and EVX2

tmp=monomers[grep("OLIG2",monomers$ID),] #These are not representative


#what represents OLIG2_EVX2

indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))

#choose the one with higher IC, "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420"
representative=metadata %>% filter(ID %in% names(represented_by_representatives)[indices]) %>%filter(IC==max(IC)) %>% pull(ID) 

#representative=names(represented_by_representatives)[indices]

tmp=metadata[which(metadata$ID %in% representative),]

monomers[grep("OLIG2",monomers$ID),] #This is not representative
#Common part of OLIG2 names

motif=find_common_prefix(monomers$ID[grep("OLIG2",monomers$ID)])
#test if any is representative
#if not then
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))
#can be many, Let's take the shortes
mono_representative=names(represented_by_representatives)[indices]
motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) #Atoh1

motif=find_common_prefix(monomers$ID[grep("EVX2_HT",monomers$ID)])
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x)))) #ALX4
mono_representative=names(represented_by_representatives)[indices]
motif2=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) 
motif2=motif2[1]



indices <- which(
  (tomtom$Query_ID==motif1) & (tomtom$Target_ID==representative) & 
    tomtom$Query_ID != tomtom$Target_ID
)

tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")] #10 9 +

indices <- which(
  (tomtom$Query_ID==motif2) & (tomtom$Target_ID==representative) & 
    tomtom$Query_ID != tomtom$Target_ID
)

tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")] #0 8 +

#In this case the coordinates are: 
#motif1 (Atoh1): start: 1+10=11 end:10(start)+9(overlap)=19
#motif2 (ALX4): start: 1+0 end: 0(start)+8 (overlap)=8
#cap position: end(ALX4)+1:start(Atoh1  )
#cap size: start(Atoh1)-end(ALX4)= 11 - 8-1=2

motif=metadata$ID[grep("FOXI1_ELF1",metadata$ID)] #not representative
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))

#choose the one with higher IC
representative2=metadata %>% filter(ID %in% names(represented_by_representatives)[indices]) %>% pull(ID) 

motif=find_common_prefix(monomers$ID[grep("FOXI1",monomers$ID)])
#test if any is representative, there is one representative dimer, do not consider that

#if not then
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))
#can be many, Let's take the shortes
mono_representative=names(represented_by_representatives)[indices]
motif21=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) #FOXA1

motif=find_common_prefix(monomers$ID[grep("ELF1",monomers$ID)]) #none of these are representative
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x)))) #
mono_representative=names(represented_by_representatives)[indices]
motif22=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) 

indices <- which(
  (tomtom$Query_ID %in% c(representative2, motif21) & tomtom$Target_ID %in% c(representative2, motif21)) & 
    tomtom$Query_ID != tomtom$Target_ID
)
#Optimal offset 3, orientation - Overlap 10

indices <- which(
  (tomtom$Query_ID %in% c(representative2, motif22) & tomtom$Target_ID %in% c(representative2, motif22)) & 
    tomtom$Query_ID != tomtom$Target_ID
)


tomtom[indices,]
#Optimal offset 4 orientation + overlap 9

#motif1(FOXA1): 13-3-10+1    13:length(dimer)-3(offset)=10 (). the match is reverse complement
#motif2(ELK1): 1+4 - 4+9

#The motifs are overlapping, this is composite

#cap size: start(ELK1)-end(FOXA1)= 5 - 10=-5
#The size of the overlapping part is abs(-5)+1=6
#Coordinates of the overlapping part: start(ELK1)-end(FOXA1)

#Mihin motif1 ja motif2 matchaa representative:ssa

#Extract the tomtom output between these motifs


sim_lt=as.data.frame(pivot_wider(tomtom, id_cols="Query_ID",
                                 names_from="Target_ID", values_from="E-value"))

sim_lt=sim_lt[, c("Query_ID", sim_lt$Query_ID)]

#sim_lt=sim_lt[,-1]
dim(sim_lt)


sim = matrix(0,dim(metadata)[1], dim(metadata)[1])
dim(sim)
tril_indices=lower.tri(sim, diag = TRUE)
ltril_indices=lower.tri(matrix(0, nrow(sim_lt),nrow(sim_lt)), diag = TRUE)


sim[tril_indices]=sim_lt[,-1][ltril_indices]
sim_sym = sim + t(sim) - diag(diag(sim), nrow(sim), ncol(sim))

dimnames(sim_sym)=list(sim_lt$Query_ID, sim_lt$Query_ID)

sim_sym[representative, motif1] #these are the same
sim_sym[motif1,representative]

sim_sym[representative, motif2] #these are the same
sim_sym[motif2,representative]



sim_df=data.frame(sim_sym, check.names = FALSE)



#tomtom ei ole välttämättä löytänyt päällekkäisyyttä monomeerin ja heterodimeerin välillä

indices <- which(
  (tomtom$Query_ID %in% c(representative, motif1) & tomtom$Target_ID %in% c(representative, motif1)) & 
    tomtom$Query_ID != tomtom$Target_ID
)

tomtom[indices,]

which(tomtom$Query_ID==representative)
which(tomtom$Target_ID==representative)





#There is two representatives, take the shorter, the other is likely multimeric
motif2=monomers[grep("FOXK1",monomers$ID)[2],] %>% pull(ID)



grep("OLIG2", capselex$ID)

grep("FOXI1_ELF2", capselex$ID)
tmp=metadata[grep("OLIG2_HT", metadata$ID),] #this is not representative


motif="ELF1_FOXK1_TTCATT40NTAT_YXI_NAGAAAACCGAAWMN_m2_c3b0_short_composite_new" #the individual motifs do not match this

#what represents ELF1_FOXK1

indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))

representative=names(represented_by_representatives)[indices]

#monomeric motifs for ELF1 and FOXK1

monomers$ID[grep("ELF1",monomers$ID)] #This is not representative
#Common part of ELF1 names




# Test the function

motif=find_common_prefix(monomers$ID[grep("ELF1",monomers$ID)])
#test if any is representative
#if not then
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))
#can be many, Let's take the shortes
mono_representative=names(represented_by_representatives)[indices]
motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID)

motif=find_common_prefix(monomers$ID[grep("FOXK1",monomers$ID)])
#There is two representatives, take the shorter, the other is likely multimeric
motif2=monomers[grep("FOXK1",monomers$ID)[2],] %>% pull(ID)






#read pwms

pcm=as.matrix( read.table(  metadata$filename[metadata$ID==representative] , header=FALSE) )

pcm1=as.matrix( read.table(  metadata$filename[metadata$ID==motif1], header=FALSE) )
pcm2=as.matrix( read.table(  metadata$filename[metadata$ID==motif2], header=FALSE) )
dimnames(pcm)=list(c("A", "C", "G", "T"))
dimnames(pcm1)=list(c("A", "C", "G", "T"))
dimnames(pcm2)=list(c("A", "C", "G", "T"))

pcm_class <- TFBSTools::PFMatrix(strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 profileMatrix=pcm
)

pcm1_class <- TFBSTools::PFMatrix(strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 profileMatrix=pcm1
)

pcm2_class <- TFBSTools::PFMatrix(strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 profileMatrix=pcm2
)


pwm_class <- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
pwm1_class <- TFBSTools::toPWM(pcm1_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
pwm2_class <- TFBSTools::toPWM(pcm2_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

#Extract the tomtom output between



source("~/projects/TFBS/RProjects/TFBS/code/KL_divergence_between_motifs.R")

compare_transforms( pwm1_class@profileMatrix,pwm_class@profileMatrix)


which(metadata$ID==mono_representative)


lapply(represented_by_representatives, function(x) grep)


#composite motifs: FIGLA-ETV1/2, HOXC10-TBX4, and FOXO1-ELK3, GLI3:RFX3
FOXK1-ELF1 composite 

TWISTTFC4/ALX3 FOXO1/ETS1, GLI3-RFX3 HOXB12/HOXD12-SOX11, HOXB2-PROX2

protein dimers (MYC:MAX, CEBPD:ATF4, FOS:JUN), DNA facilitated dimers (MEIS1:HOXB13, MEIS1:DLX3, ETS1 homodimer) 
and DNA dependent dimers (ETS1:PAX5, ERG:FLI1?, some other) 17. 
The structures included BARHL2 dimer with two spacings, TEAD4 dimer, HOXB13 dimer, 
TEAD4-HOXB13 and MEIS1-HOXB13 heterodimers, and FOXK1-ELF1 heterodimer.

TEAD4-HOXB13



artificial_all<-list()

for(i in 1:nrow(metadata) ){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
  
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  min_slice_width = floor( ncol(pcm) / 3 )
  
  artificial<-list()
  for( j in min_slice_width:(ncol(pcm) - min_slice_width) ){
    #j=5
    # Get matrix versions where slices are concatenated in all 3 relative orientations
   
    artificial[paste0(metadata$ID[i], "_",j,c("_HT2", "_HH", "_TT")) ]=get_other_pfm_slice_pair_orientations(pcm, j)
    
  }
  
  for(motif in names(artificial)){
    print(motif)
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../artificial_motifs/", motif, ".pfm"), sep="\t")
    write.table(artificial[[motif]],row.names = FALSE, col.names=FALSE, quote=FALSE,file=paste0("../../artificial_motifs_space/", motif, ".pfm"), sep=" ")
  }
  
  artificial_all[[metadata$ID[i]]]=artificial
  
  #Write all motifs to file:
  

}

saveRDS(artificial_all, file="Rdata/artificial_all.Rds")

stop()



data_path="/scratch/project_2006203/"

#df_motif_info <- read_tsv("../metadata/df_motif_info.tsv", col_names=TRUE)
#representatives <- read_tsv("../metadata/new_representatives.tsv", col_names=TRUE)
representatives <- read_tsv(paste0(data_path,"motif-clustering-Viestra-private/metadata/new_representatives_IC_length_Morgunova.tsv"), col_names=TRUE)

for(i in 1:nrow(representatives)){
  #i=1
  print(i)
  pcm=as.matrix( read.table(paste0(data_path,"TFBS/", representatives$filename[i]), header=FALSE) )
  length=ncol(pcm)
  
  #In how many different combinations you can have #length columns
  factorial(length) #40320
  
  #Generate all permutations
  library(gtools)
  perms <- gtools::permutations(ncol(pcm), ncol(pcm))
  perms <- gtools::permutations(ncol(pcm), ncol(pcm))
  pairings <- expand.grid(1:col(perms), 1:ncol(perms))
  
  all.perms <- lapply(1:nrow(pairings), function(x) pcm[, perms[pairings[x,1],] ])
  
  all.unique.perms <- unique(perms)
  length(all.unique.perms)
  
  perms2=Permn(1:length)
  perms=perms[-1,] #the original not wanted
  
  motif_perms=pcm[,perms[1,]]
  
  #For creating the whole dataset of the permutations, there’s the function DescTools::CombSet(), which has
  #an argument m for defining the size of the subset to be drawn and where the replacement and order
  #arguments can be set. 
  #CombN(1:length, 4, repl=FALSE, ord=FALSE) 

  #Sample two or three columns at a time
  
  #How many two column parts?
  
  #length(which((representatives$length %% 2) == 0)) #even 4218
  #length(which((representatives$length %% 2) != 0)) #odd 1565
  
  if(length %% 2==0){ #even
    
    parts=length/2 #length is even
    
    ind=as.matrix(data.frame(first=seq(1,(length-1),2), second=seq(2,length,2)))
    colnames(ind)=NULL
    ind.list <- split(ind, seq(nrow(ind)))
    
    ind2=as.matrix(data.frame(first=seq(2,length,2), second=c( seq(3,length-1,2),1) ))
    colnames(ind2)=NULL
    ind2.list <- split(ind2, seq(nrow(ind2)))
    
  }else{ #odd
    
    parts=ceiling(length/2)
    
    #first position alone
    ind=as.matrix(data.frame(first=c(1, seq(2,(length-1),2)), second=c(NA,seq(3,length,2))) )
    colnames(ind)=NULL
    ind.list <- split(ind, seq(nrow(ind)))
    ind.list[[1]]=1
    
    #last position alone
    ind2=as.matrix(data.frame(first=c(seq(1,length-2,2),length), second=c( seq(2,length-1,2),length) ))
    colnames(ind2)=NULL
    ind2.list <- split(ind2, seq(nrow(ind2)))
    ind2.list[[parts]]=length
    
  }
  
  #How many combinations of parts
  perms=Permn(1:parts)
  perms=perms[-1,] #the original not wanted
  
  #all possible motif permutations
  
  

  
  #reverse complements
  #A control set of reversed but not complemented motifs was generated by reversing the column order of each motif matches (Sahu et al. 2022)
  
  #handle cap-selex motifs
  

}

# x <- letters[1:4]
# n <- length(x)
# m <- 2
# factorial(n) #24
# Permn(x) 
# CombN(n, m, repl=FALSE, ord=TRUE) #12
# CombSet(x, m, repl=FALSE, ord=TRUE) #All combinations of two letters, the letter can not occur twice, the order matters
# 
# CombN(n, m, repl=TRUE, ord=TRUE) #16
# CombSet(x, m, repl=TRUE, ord=TRUE) #All combinations of two letters, the letter can not occur twice, the order matters
# 
# CombN(n, m, repl=TRUE, ord=FALSE) #10
# CombSet(x, m, repl=TRUE, ord=FALSE) #All combinations of two letters, the letter can occur twice, the order does not matter
# 
# CombN(n, m, repl=FALSE, ord=FALSE) #6
# CombSet(x, m, repl=FALSE, ord=FALSE) #All combinations of two letters, the letter can not occur twice, the order does not matter



df_motifsimilarity=read_tsv('../metadata/motifsimilarity.tsv') #7_926_171 rows
#tmp=df_motifsimilarity.tail(1829)

n=3982
k=2
P = choose(n, k)
  
factorial(n) / ( factorial(k)* factorial(n - k) )
print(P) #n*(n-1)/2=7_926_171


min(df_motifsimilarity$gapped)
max(df_motifsimilarity$gapped)

library(patchwork)
#install.packages("patchwork")
library(ggbreak)
#install.packages("ggbreak")
pdf("Figures/histogram-of-gapped-similarities.pdf", width=4, height=4)
ggplot(df_motifsimilarity, aes(x=gapped)) +  geom_histogram(aes(y = after_stat(count / sum(count))), bins = 200) +
  scale_y_continuous(labels = scales::percent)+ scale_y_cut(breaks=c(0.000075, 0.01, 0.93), which=c(1,2,3,4), scales=c(0.5,0.5,0.5,3))+
  #scale_y_break(c(0.0225, 0.93 ))
  #geom_histogram(color="black", fill="white", bins=100)+xlim(0,1)+
  labs(x='Similarity', y='Frequency', title='Histogram of gapped 10-mer similarities')
dev.off()



sim_lt=pivot_wider(df_motifsimilarity, id_cols="Query_ID",
                   names_from="Target_ID", values_from="gapped")

half_lt=pivot_wider(df_motifsimilarity, id_cols="Query_ID",
                   names_from="Target_ID", values_from="half")


sim = matrix(0,dim(representatives)[1], dim(representatives)[1])
half= matrix(0,dim(representatives)[1], dim(representatives)[1])

dim(sim)
dim(half)
diag(sim)=1
diag(half)=1

tril_indices=lower.tri(sim, diag = FALSE)

ltril_indices=lower.tri(matrix(0, nrow(sim_lt),nrow(sim_lt)), diag = TRUE)

#sim[row, col]=df_motifsimilarity.gapped
sim[tril_indices]=sim_lt[,-1][ltril_indices]
half[tril_indices]=half_lt[,-1][ltril_indices]

sim_sym = sim + t(sim) - diag(diag(sim), nrow(sim), ncol(sim))
half_sym = half + t(half) - diag(diag(half), nrow(half), ncol(half))

a=colnames(sim_lt)[2]
  

b=sim_lt[,"Query_ID"]$Query_ID

a=c(a,b)

rownames(sim_sym)=a
colnames(sim_sym)=a

rownames(half_sym)=a
colnames(half_sym)=a


sim_df=data.frame(sim_sym, check.names = FALSE)
dist_df =data.frame(1-sim_sym, check.names=FALSE)

half_df=data.frame(half_sym, check.names = FALSE)
half_dist_df =data.frame(1-half_sym, check.names=FALSE)


#Some motif names in sim_df and representatives do not match?

#missing=which(!(names(sim_df) %in% representatives$ID))
#missing=which(!( representatives$ID %in% names(sim_df)))
#newnames=names(sim_df)
#newnames[missing]=paste0(newnames[missing] ,"_NA")


#rownames(sim_df)=newnames
#colnames(sim_df)=newnames

#rownames(dist_df)=newnames
#colnames(dist_df)=newnames




#Consider 1062 representatives

table(representatives$new_representative)

rep_motifs=representatives %>% filter(new_representative=="YES") %>% select(ID)

sim_rep_df=sim_df[rep_motifs$ID, rep_motifs$ID]
dist_rep_df=1-sim_rep_df

half_rep_df=half_df[rep_motifs$ID, rep_motifs$ID]
half_dist_rep_df=1-half_rep_df


#Cluster the matrix, consider rows as samples, columns as features 



#CAP-selex motifs for which both TFs from the same protein family

CAP_index=which(representatives$experiment=="CAP-SELEX")

tmp=representatives[CAP_index, "Lambert2018.families"]

res = strsplit(tmp[[1]], "_")
str(res)  

res=do.call(rbind, res)

CAP_index=CAP_index[which( res[,1]==res[,2])] #120

representatives <- representatives %>% mutate(family_alt=`Lambert2018.families`)

representatives[CAP_index, "family_alt"]=res[ res[,1]==res[,2],1]

tmp=representatives %>% count(family_alt, sort=TRUE)



rm=unique( c( grep("_", tmp[[1]]), grep(";", tmp[[1]])) )

tmp=tmp[-rm,]

#choose only first 10

tmp=tmp[1:10,"family_alt"]

#These need to be in same order as in sim_df
anno=representatives %>% select(ID, "family_alt")  #N times D 
#remove those family_alts that do not match tmp
length(which( !(anno$family_alt %in% tmp$family_alt))) #1684
anno$family_alt[!(anno$family_alt %in% tmp$family_alt)]=NA
anno=anno[order(match(anno$ID,rownames(sim_df))),]

#are the rownames of annos and rownames of sim_df the same

mat=sim_df
min(mat)
max(mat)

#correlation does not work here
#ed <- dist(sim_df,method="correlation") #Euclidean is the default, N x D matrix as input

# Pairwise correlation between samples (columns)
cols.cor <- cor(sim_df, use = "pairwise.complete.obs", method = "pearson")
cols.dist=as.dist(1 - cols.cor)

hc_cor_complete <- cluster::agnes(cols.dist, diss=TRUE, method="complete")

#hc<-hclust(as.matrix(sim_df), method="complete") does not work with big data
hc_complete <- cluster::agnes(as.matrix(sim_df), diss=TRUE, method="complete")

plot(as.dendrogram(hc_cor_complete), which.plots=2)


fit<-kmeans(df_top1500,4, nstart=25)
fviz_cluster(fit, df_top1500,ellipse.type="norm") +theme_minimal()
#library(cluster)
sil1<-silhouette(fit$cluster,dist(df_top1500))

plot(silhouette(cutree(hc_cor_complete,4),cols.dist,border = NA))

silhouette_score <- function(k){
  print(k)
  ss<-silhouette(cutree(hc_cor_complete,k), cols.dist,border = NA)
  mean(ss[,3])
  
}
k=2:1500
k=850:1000
agv_sil <- sapply(k, silhouette_score) #MAX IS 914
plot(k, type="b", agv_sil, xlab="Number of clusters", ylab="Average silhouette score", frame=FALSE)

plot(silhouette(cutree(cluster_complete,5), cluster_distance), col = c("blue","red","purple","green","black"), border =NA)


res<-NbClust(data=mat, diss = cols.dist, distance=NULL,min.nc=2, max.nc=5, 
             method = "complete", index = "all")

data<-iris[,-c(5)]
diss_matrix<- dist(data, method = "euclidean", diag=FALSE)
diss_matrix<- factoextra::get_dist(sim_df, method = "pearson", diag=FALSE)
diss_matrix<- factoextra::get_dist(sim_df, method = "pearson", diag=TRUE)


clust.obj=NbClust(data, diss=diss_matrix, distance = NULL, min.nc=2, max.nc=6, 
        method = "complete", index = "alllong")  

#hcut is from factoextra, you can give it hc_func="hclust/agnes/diana"
fviz_nbclust(sim_df, FUN=hcut, hc_func="agnes", hc_method="complete", hc_metric="pearson", method = "silhouette", min.nc=800, max.nc=1200)+ theme_classic()
fviz_nbclust(sim_df, FUN=hcut, hc_func="agnes", hc_method="complete", hc_metric="pearson", method = "silhouette", k.max=5)+ theme_classic()

wss.plot=fviz_nbclust(sim_df, FUN=hcut, diss=diss_matrix,hc_method="complete", hc_metric=NULL, method = "wss", k.max=5)+ theme_classic()

silhouette.plot=fviz_nbclust(sim_df, FUN=hcut, diss=diss_matrix,hc_method="complete", hc_metric="pearson", method = "silhouette", k.max=20)+ theme_classic()

gap_stat <- clusGap(dsim_df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

for_plot <-  sim_df %>%
  as.data.frame() %>%
  dist(method = "correlation") %>%
  hclust(method = "complete") %>%
  cutree(h=25, k=2)

sil<- for_plot %>%
  silhouette(.,distance_mat)
plot(sil)

#dist.euc <- daisy(dat, metric = "euclidean")
#single <- agnes(dist.euc, method = "single")
#plot(as.dendrogram(single), leaflab = "none",
#     main = "Single linkage dendrogram")


#library('factoextra')
subset_mat <- t(subset_mat)
distance_mat <- dist(subset_mat, method = 'euclidian')
fviz_nbclust(subset_mat, FUNcluster = hcut ,method = "silhouette")





png("Figures/motifsimilarity-heatmap_all_complete_correlation.png", width=3100,height=2000, res=315, pointsize=12)
ComplexHeatmap::Heatmap(mat, 
         name="Gapped 10-mer similarity",               
         clustering_distance_rows="pearson",
         clustering_distance_columns="pearson", 
         col= circlize::colorRamp2( c(0, 1), c( "white", "blue") ),
         clustering_method_rows="complete",
         clustering_method_columns="complete",
         #cluster_rows=rev(as.dendrogram(hc)),
         #cluster_columns=rev(as.dendrogram(hc)),
         show_row_names=FALSE,
         show_column_names=FALSE, 
         left_annotation = rowAnnotation(family=anno$family_alt, 
                                             col=list(family=c("Homeodomain"="blue", "C2H2 ZF"="#FF7F00", "bHLH"="#66A61E",
                                                               "bZIP"="#E41A1C", "Nuclear receptor"="purple", "Ets"="brown", 
                                                               "HMG/Sox"="#E7298A", "Forkhead"="gray", "T-box"="#FFE528", "IRF"="cyan", "NA"="white")
                                               )), 
         show_row_dend=TRUE, show_column_dend = FALSE,
         row_dend_width=unit(50, "pt"), 
         #column_dend_height=unit(0,"pt"),
         height = nrow(mat)*unit(0.1,"pt"), 
         width = ncol(mat)*unit(0.1, "pt"), 
         use_raster=TRUE,
          #file="Figures/motifsimilarity-heatmap_all_complete_correlation.png"
)
dev.off()

#half and gapped in the same figure


sim_both = matrix(0,dim(sim_df)[1], dim(sim_df)[1])
tril_indices=lower.tri(sim_both, diag = TRUE)
triu_indices=upper.tri(sim_both, diag = FALSE)

sim_both[tril_indices]=sim_df[tril_indices]
sim_both[triu_indices]=half_df[triu_indices]

utril_indices=upper.tri(matrix(0, nrow(sim_df),nrow(sim_df)), diag = FALSE)

r = ((triu_indices-1) %% nrow(sim_df)) + 1
c = floor((triu_indices-1) / nrow(sim_df)) + 1


mat[triu_indices]=as.matrix(half_df)[triu_indices]

mat[r,c]=half_df[r,c]

mat=as.matrix(sim_both)
mat=mat[rev(hc$order), rev(hc$order)]

test=mat[1:5,1:5]
rownames(test)=NULL
colnames(test)=NULL

zeromat=matrix(0, nrow=5, ncol=5)
triu_indices=upper.tri(as.matrix(test), diag = FALSE)
zeromat[triu_indices]=test[triu_indices]

png("Figures/gapped-half-heatmap_all_ward.png", width=3100,height=2000, res=315, pointsize=12)
ComplexHeatmap::Heatmap(mat, 
                        name="Gapped 10-mer similarity",               
                        #clustering_distance_rows="pearson",
                        #clustering_distance_columns="pearson", 
                        col= circlize::colorRamp2( c(0, 1), c( "white", "blue") ),
                        #clustering_method_rows="complete",
                        #clustering_method_columns="complete",
                        cluster_rows=FALSE, #as.dendrogram(hc),
                        cluster_columns=FALSE, #as.dendrogram(hc),
                        #cluster_rows=FALSE,
                        #cluster_columns=FALSE,
                        show_row_names=FALSE,
                        show_column_names=FALSE, 
                        left_annotation = rowAnnotation(family=anno$family_alt[rev(hc$order)], 
                                                        col=list(family=c("Homeodomain"="blue", "C2H2 ZF"="#FF7F00", "bHLH"="#66A61E",
                                                                          "bZIP"="#E41A1C", "Nuclear receptor"="purple", "Ets"="brown", 
                                                                          "HMG/Sox"="#E7298A", "Forkhead"="gray", "T-box"="#FFE528", "IRF"="cyan", "NA"="white")
                                                        )), 
                        show_row_dend=FALSE, show_column_dend = FALSE,
                        #row_dend_width=unit(50, "pt"), 
                        #column_dend_height=unit(0,"pt"),
                        height = nrow(mat)*unit(0.1,"pt"), 
                        width = ncol(mat)*unit(0.1, "pt"), 
                        use_raster=TRUE,
                        #file="Figures/motifsimilarity-heatmap_all_complete_correlation.png"
)
dev.off()

mat=as.matrix(sim_both)
mat=mat[hc_complete$order, hc_complete$order]

png("Figures/gapped-half-heatmap_all_complete.png", width=3100,height=2000, res=315, pointsize=12)
ComplexHeatmap::Heatmap(mat, 
                        name="Half and gapped 10-mer similarity",               
                        #clustering_distance_rows="pearson",
                        #clustering_distance_columns="pearson", 
                        col = viridisLite::viridis(n = 256,  option = "magma"),
                        #col= circlize::colorRamp2( c(0, 1), c( "white", "blue") ),
                        #clustering_method_rows="complete",
                        #clustering_method_columns="complete",
                        cluster_rows=FALSE, #as.dendrogram(hc),
                        cluster_columns=FALSE, #as.dendrogram(hc),
                        #cluster_rows=FALSE,
                        #cluster_columns=FALSE,
                        show_row_names=FALSE,
                        show_column_names=FALSE, 
                        left_annotation = rowAnnotation(family=anno$family_alt[hc_complete$order], 
                                                        col=list(family=c("Homeodomain"="blue", "C2H2 ZF"="#FF7F00", "bHLH"="#66A61E",
                                                                          "bZIP"="#E41A1C", "Nuclear receptor"="purple", "Ets"="brown", 
                                                                          "HMG/Sox"="#E7298A", "Forkhead"="gray", "T-box"="#FFE528", "IRF"="cyan", "NA"="white")
                                                        ), annotation_legend_param=list(at=tmp$family_alt)), 
                        show_row_dend=FALSE, show_column_dend = FALSE,
                        #row_dend_width=unit(50, "pt"), 
                        #column_dend_height=unit(0,"pt"),
                        height = nrow(mat)*unit(0.1,"pt"), 
                        width = ncol(mat)*unit(0.1, "pt"), 
                        use_raster=TRUE,
                        #file="Figures/motifsimilarity-heatmap_all_complete_correlation.png"
)
dev.off()


png("Figures/gapped-heatmap_all_complete.png", width=3100,height=2000, res=315, pointsize=12)
ComplexHeatmap::Heatmap(as.matrix(sim_df), 
                        name="Gapped 10-mer similarity",               
                        #clustering_distance_rows="pearson",
                        #clustering_distance_columns="pearson", 
                        col= circlize::colorRamp2( c(0, 1), c( "white", "blue") ),
                        #clustering_method_rows="complete",
                        #clustering_method_columns="complete",
                        cluster_rows=as.dendrogram(hc_complete),
                        cluster_columns=as.dendrogram(hc_complete),
                        #cluster_rows=FALSE,
                        #cluster_columns=FALSE,
                        
                        show_row_names=FALSE,
                        show_column_names=FALSE, 
                        left_annotation = rowAnnotation(family=anno$family_alt, 
                                                        col=list(family=c("Homeodomain"="blue", "C2H2 ZF"="#FF7F00", "bHLH"="#66A61E",
                                                                          "bZIP"="#E41A1C", "Nuclear receptor"="purple", "Ets"="brown", 
                                                                          "HMG/Sox"="#E7298A", "Forkhead"="gray", "T-box"="#FFE528", "IRF"="cyan", "NA"="white")
                                                        ),
                                                        annotation_legend_param=list(at=tmp$family_alt)
                                                        ), 
                        show_row_dend=TRUE, show_column_dend = FALSE,
                        row_dend_width=unit(100, "pt"), 
                        #column_dend_height=unit(0,"pt"),
                        height = nrow(mat)*unit(0.1,"pt"), 
                        width = ncol(mat)*unit(0.1, "pt"), 
                        use_raster=TRUE,
                        #file="Figures/motifsimilarity-heatmap_all_complete_correlation.png"
)
dev.off()





# png("Figures/motifsimilarity-heatmap_all_complete_correlation.png", width=3100,height=2000, res=315, pointsize=12)
pheatmap(as.matrix(sim_df), clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", color=color,clustering_method="complete",
         show_rownames=FALSE,
         show_colnames=FALSE, annotation_row=annos, treeheight_row=0, treeheight_col=0,
         cellheight = 0.1, cellwidth = 0.1, 
         annotation_colors=ann_colors #, #cellheight = 1, cellwidth = 1,
         #file="Figures/motifsimilarity-heatmap_all_complete_correlation.png"
         )
dev.off()



#pdf("Figures/heatmap_all_ward_correlation.pdf")
png("Figures/motifsimilarity-heatmap_all_ward_correlation.png", width=3100,height=2000, res=315, pointsize=12)
pheatmap(sim_df, clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", color=color,clustering_method="ward.D",
         show_rownames=FALSE, treeheight_row=0, treeheight_col=0,
         cellheight = 0.1, cellwidth = 0.1, 
         show_colnames=FALSE, annotation_row=annos, annotation_colors=annoCol #, 
         #file="Figures/motifsimilarity-heatmap_all_ward_correlation.png") 
)
dev.off()

# annotation=anno, clustering_method="complete",
#         drop_levels=TRUE, legend=FALSE,
#         cutree_cols=2)



ehc <- hclust(ed,method="complete")
ehc2 <- hclust(ed,method="single")
# plot the resulting dendrogram: complete linkage method produces clear groups (2 or 3 clear clusters), while single linkage starts from a single patient and goes to the most dissimilar which is not very informative here.
plot(ehc,labels=FALSE)
plot(ehc2, labels=F)


