

rm(list=ls())

source("code/half_site_recognition_functions.R", echo=FALSE)

library(readr)
library(tidyverse)

#TFBS_path='/Users/osmalama/projects/TFBS/'


tomtom = read_tsv('/Users/osmalama/projects/TFBS/Results/tomtom_final_version2.2/tomtom_relaxed/tomtom.tsv',
                   trim_ws = TRUE) # 15 468 492 rows

#remove three last rows, this was the source of parsing issues reported by the previos command
tomtom <- tomtom %>% dplyr::slice(1:(n()-3)) #15 468 489

# Distances between artificial HT-SELEX half-site motifs and version2.2. motifs

tomtom_halfsites = read_tsv('/Users/osmalama/projects/TFBS/Results/tomtom_final_version2.2/tomtom_relaxed_artificialHTSelex_against_version2.2/tomtom.tsv') # 15 468 492 rows

#remove three last rows
tomtom_halfsites <- tomtom_halfsites %>% dplyr::slice(1:(n()-3)) #27 531 = 7 * 3933 (Correct!)


tomtom=rbind(tomtom, tomtom_halfsites) #154 960 20 

#In reality 3933 motifs, typo in filename
metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

#Metadata for artificial half sites, 7 motifs
metadata_halfsites <- read_delim("~/projects/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)


#Which monomers need to be split/get monomeric motif
# comments_from_curation<- read_csv("~/projects/TFBS/PWMs_final_union/metadata_google_drive_supplement_manual_curation_spacing_composite.csv")
# unique(comments_from_curation$IS_commets)


# tmp=do.call(rbind, sapply(as.list(comments_from_curation$symbol[unique(c(grep(c("homodimer"),comments_from_curation$IS_commets),
#   grep(c("Homodimer"),comments_from_curation$IS_commets),
#   grep(c("Trimer"),comments_from_curation$IS_commets),
#   grep(c("trimer"),comments_from_curation$IS_commets)
# ))]), strsplit,"_"))

# tmp=unique(c(tmp[,1], tmp[,2])) #potential homodimers

# unique(comments_from_curation$YY_comments)
# unique(comments_from_curation$JT_comments)

# tmp2=do.call(rbind, sapply(as.list(comments_from_curation$symbol[grep("halfsite",comments_from_curation$JT_comments)]), strsplit, "_"))
# tmp2=unique(c(tmp2[,1], tmp2[,2])) #potential homodimers
# 
# tmp=unique(c(tmp, tmp2))

# unique(comments_from_curation$YC_comments)
# tmp3=do.call(rbind, sapply(as.list(comments_from_curation$symbol[unique(c(grep("half-site",comments_from_curation$YC_comments),
#                               grep("homodimer",comments_from_curation$YC_comments)))]), strsplit, "_"))
# 
# tmp3=unique(c(tmp3[,1], tmp3[,2])) #potential homodimers
# tmp=unique(c(tmp, tmp3))
# 

# "MAFG" ?
# "NRL" ?
# "POU3F2" ?
# "POU4F1" ?
# "MAFA" ?
# "MAF" ?
# "PAX7" ?
# "PAX2" ?
# "POU4F3" ?
#   
# "OTP" no monomer
#   
# "IRF8" homodimer
# "IRX5" homodimer
# "MAFF" homodimer
# "NFATC1" homodimer
# "PROP1" homodimer?
# "FOXD2" homodimer
# "SOX11" homodimer
# "TGIF2LX" homodimer
# "IRF7"  homodimer
# "NFKB2" homodimer
# "RFX3" homodimer
# "ONECUT1" homodimer
# "RFX5" homodimer
# "PROX1" homodimer?
# TFAP2C: no monomer?


#monomers
monomers=metadata[metadata$experiment!="CAP-SELEX",]

monomers_halfsites=metadata_halfsites #"FOXD3"(these should be this)   "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"





#monomers_representative=monomers[which(monomers$new_representatives=="YES"),]

#How to check which motif is spacing and which is composite?

#representatives=metadata[which(metadata$new_representative=="YES"),] #metadata for 1031 representatives

#represented_by_representatives <- read_table("~/projects/TFBS/PWMs_final/sstat_represented_by_representatives.tsv")

#represented_by_representatives <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/represented_by_representatives.RDS")

#which are spacing or composite motifs
#capselex=metadata[metadata$experiment=="CAP-SELEX",]

# capselex=read_delim("~/projects/TFBS/PWMs_final_union/metadata_google_drive_supplement.tsv", 
#                     delim = "\t", escape_double = FALSE, 
#                     trim_ws = TRUE)

capselex= metadata %>% filter(experiment=="CAP-SELEX") #1898

capselex %>% filter(study=="fromYimeng") %>% select(type) %>% table()

#capselex$type[grep("composite",capselex$ID)]="composite"
#capselex$type[grep("spacing",capselex$ID)]="spacing"
#capselex$type[grep("20230420",capselex$type)]="spacing"

capselex$first.monomer=""
capselex$second.monomer=""

capselex$First_dimeric_bool=NA
capselex$First_tomtom_reverse_order=NA
capselex$First_Optimal_offset=NA
capselex$First_Optimal_Overlap=NA
capselex$First_Optimal_Orientation=""

capselex$Second_dimeric_bool=NA
capselex$Second_tomtom_reverse_order=NA
capselex$Second_Optimal_offset=NA
capselex$Second_Optimal_Overlap=NA
capselex$Second_Optimal_Orientation=""



##"FOXD3"   "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"
artificial_halfsites=list()
artificial_halfsites[["IRF8"]]="IRF6"
artificial_halfsites[["MAFF"]]="MAFF"
artificial_halfsites[["NFATC1"]]="NFATC1"
#artificial_halfsites[["FOXD2"]]="FOXD3"
artificial_halfsites[["SOX11"]]="SOX11"
artificial_halfsites[["TGIF2LX"]]="TGIF2LX"
artificial_halfsites[["ONECUT1"]]="ONECUT1"

# alternatives=list()
# alternatives[["BHLHB8"]]="BHLHA15"
# alternatives[["PROX2"]]="PROX1"
# alternatives[["NKX2-4"]]="NKX2-3"
# alternatives[["DBX2"]]="DLX2"
# alternatives[["BHLHB5"]]="BHLHA15"
# alternatives[["ZBTB33"]]="ZBTB33" #?
# alternatives[["PBX2"]]="PBX2" #?
# alternatives[["FOXI2"]]="FOXI1"
# alternatives[["gntR"]]="gntR" #?
# alternatives[["NKX2-6"]]="NKX2-3"

#tomtom fails
#808  814  828  830  831  832  836  845  892  895  945  951  980  983  998 1056 1444 1445 1575 1653 1760 1761 1762 1763 1768 1769 1770 1778 1779


#which(capselex$ID=="TCF4_FOXD2_TAAGAT40NATG_YWI_NCACCTGWAATAACAN_m1_c3b0u") 1595
missing_monomers=c()
tomtom_issue=c()

for(i in 1:nrow(capselex)){
  #i=1
  print(i)
  dimer_ID=capselex$ID[i] #Need to remove possible "-" from the ID as the tomtom names do not have these
  
  #Check whether the last character is "-" and remove it
  if(substr(dimer_ID, nchar(dimer_ID), nchar(dimer_ID))=="-"){
    dimer_ID=substr(dimer_ID, 1, nchar(dimer_ID)-1)
  }
  
  #monomers
  monomers_symbols=unlist(strsplit(capselex$symbol[i],"_"))
  
  # Find representative monomers, experimental
  # motif1=find_representative_monomer(monomer=monomers_symbols[1], monomers, represented_by_representatives, metadata)
  # motif2=find_representative_monomer(monomer=monomers_symbols[2], monomers, represented_by_representatives, metadata)
  
  #Find monomers by symbol only, for which motifs this fails.
  if(monomers_symbols[1] %in% names(artificial_halfsites)){
    motif1=list()
    motif1$motif=metadata_halfsites %>% filter(symbol==as.character(artificial_halfsites[monomers_symbols[1]])) %>% select(ID) %>% pull(ID)
    motif1$representative_bool=FALSE
    motif1$dimeric_bool=FALSE
    
    
  }else{
    motif1=try(find_monomer(monomer=monomers_symbols[1], monomers), silent=TRUE)
  }
  
  if(monomers_symbols[2] %in% names(artificial_halfsites)){
    motif2=list()
    motif2$motif=metadata_halfsites %>% filter(symbol==as.character(artificial_halfsites[monomers_symbols[2]])) %>% select(ID) %>% pull(ID)
    motif2$representative_bool=FALSE
    motif2$dimeric_bool=FALSE
  }else{
    motif2=try(find_monomer(monomer=monomers_symbols[2], monomers), silent=TRUE)
  }
  
  indices=NULL
  
  if(length(motif1$motif)==1){
  
    tomtom_reverse_order=FALSE
    indices <- which( #92250
      (tomtom$Query_ID==motif1$motif) & (tomtom$Target_ID==dimer_ID) & 
        tomtom$Query_ID != tomtom$Target_ID
    )
    
    if(length(indices)==0){
      tomtom_reverse_order=TRUE
      indices <- which( #1309620
        (tomtom$Target_ID==motif1$motif) & (tomtom$Query_ID==dimer_ID) & 
          tomtom$Query_ID != tomtom$Target_ID
      )
    }
    
    capselex$first.monomer[i]=motif1$motif
    capselex$First_dimeric_bool[i]=motif1$dimeric_bool
    capselex$First_tomtom_reverse_order[i]=tomtom_reverse_order
    
    
    if(length(indices)!=0){
      capselex[i, c("First_Optimal_offset","First_Optimal_Overlap","First_Optimal_Orientation")]=tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")]
    }else{
      print("No tomtom match found for monomer1")
      tomtom_issue=c(tomtom_issue, motif1$motif)
    }
  
  }else{
    print("No first monomer found")
    missing_monomers=c(missing_monomers, monomers_symbols[1])
  }
  
  if(length(motif2$motif)==1){
  
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
  
  
  capselex$second.monomer[i]=motif2$motif
  capselex$Second_dimeric_bool[i]=motif2$dimeric_bool
  capselex$Second_tomtom_reverse_order[i]=tomtom_reverse_order
  
  if(length(indices)!=0){
    capselex[i, c("Second_Optimal_offset","Second_Optimal_Overlap","Second_Optimal_Orientation")]=tomtom[indices,c("Optimal_offset", "Overlap", "Orientation")]
  }else{
    print("No tomtom match found for monomer1")
    tomtom_issue=c(tomtom_issue, motif2$motif)
  }
  
  
  }else{
    print("No second monomer found")
    missing_monomers=c(missing_monomers, monomers_symbols[2])
  }
  
}


print(unique(missing_monomers))
#"PBX4"   "BHLHB8" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" "BHLHB5" "NKX2-6"
print(unique(tomtom_issue))

save.image(file="~/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs.RData")
load(file="~/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs.RData")

#"BHLHB8" "PBX4" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" "BHLHB5" "NKX2-6"

capselex %>% filter(type=="composite")

#For how many CAP-selex motifs, either TF of the pair are in missing_monomers
pairs=as.data.frame(do.call(rbind, strsplit(capselex$symbol, split="_")))
names(pairs)=c("TF1", "TF2")
length(which(pairs$TF1 %in% unique(missing_monomers) | pairs$TF2 %in% unique(missing_monomers))) #31, 1.6% of all CAP-selex motifs

tmp=capselex[(which(pairs$TF1 %in% unique(missing_monomers) | pairs$TF2 %in% unique(missing_monomers))),] #31

table(capselex$first.monomer=="" | capselex$second.monomer=="" )
#FALSE  TRUE 
#1867    31

#For how many CAP-selex motifs, there are individual TF motifs but tomtom fails, NONE for relaxed set 

table(capselex$first.monomer!="" & capselex$second.monomer!="" )
#FALSE  TRUE 
#31  1864 

#There are 1867 for which both monomers were found
#For how many of these 1245 motifs, tomtom did not find the match for either of the monomers  
#For the relaxed setting, all monomers are found
table( (capselex$first.monomer!="" & capselex$second.monomer!="" ) & (is.na(capselex$First_Optimal_offset) | is.na(capselex$Second_Optimal_offset)))
#FALSE 
#1898 



#In total there are 26+202=389? motifs for which there are no full info for the alignment
table(is.na(capselex$First_Optimal_offset) | is.na(capselex$Second_Optimal_offset))
#FALSE  TRUE 
#1867    31 

#Extract shared columns

shared_columns <- intersect(names(metadata), names(metadata_halfsites))

#Combine the two data.frames

metadata <- rbind(metadata[shared_columns], metadata_halfsites[, shared_columns])




#match for both motifs
capselex$first_monomer_start=NA
capselex$first_monomer_end=NA

capselex$second_monomer_start=NA
capselex$second_monomer_end=NA


capselex$first_start=NA
capselex$first_end=NA

capselex$second_start=NA
capselex$second_end=NA

capselex$type_computation=""
capselex$overlap_length=NA
capselex$overlap_start=NA
capselex$overlap_end=NA

capselex$gap_length=NA
capselex$gap_start=NA
capselex$gap_end=NA

#Figure out also the position of the overlap in both monomers

which(capselex$ID=="ATF3_ONECUT1_TCATTC40NATT_YJII_NTGACGTCANNNNATCGATN_m1_c3b0")
which(capselex$ID=="MAFA_FOSB_TTACTT40NTAT_YWIII_YGCTGASTCAY_m1_c3b0")

for (i in which( paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation) %in% c("- -", "- +", "+ +", "+ -")) ){
  print(i)
  #i=which( paste0(capselex$First_Optimal_Orientation, " ", capselex$Second_Optimal_Orientation) %in% c("- -", "- +", "+ +", "+ -"))[1]
  #"first.monomer"               "second.monomer"              "First_representative_bool"   "First_tomtom_reverse_order"  "First_Optimal_offset"       
  #"First_Optimal_Overlap"       "First_Optimal_Orientation"   "Second_representative_bool"  "Second_tomtom_reverse_order" "Second_Optimal_offset"      
  #"Second_Optimal_Overlap"      "Second_Optimal_Orientation" 
  
  #which(capselex$ID =="ATOH1_ONECUT1_TCCATG40NCAC_YJII_CATATGNNNNNNATCGAT_m1_c3b0_short_spacing_new")
  
  
  first_se=monomer_start_end( tomtom_reverse_order=capselex$First_tomtom_reverse_order[i], #FALSE
                     optimal_offset=capselex$First_Optimal_offset[i], #12
                     optimal_orientation=capselex$First_Optimal_Orientation[i], #-
                     optimal_overlap=capselex$First_Optimal_Overlap[i], #8
                     dimer_length=capselex$length[i]) #25
  
  
  second_se=monomer_start_end( tomtom_reverse_order=capselex$Second_tomtom_reverse_order[i], #FALSE
                              optimal_offset=capselex$Second_Optimal_offset[i], #-1
                              optimal_orientation=capselex$Second_Optimal_Orientation[i], #+#
                              optimal_overlap=capselex$Second_Optimal_Overlap[i], #9
                              dimer_length=capselex$length[i])
  
  
  reverse_bool=list()
  reverse_bool[["TRUE"]]=FALSE
  reverse_bool[["FALSE"]]=TRUE
  
  first_monomer_se=monomer_start_end( tomtom_reverse_order=reverse_bool[[as.character(capselex$First_tomtom_reverse_order[i])]], #FALSE
                              optimal_offset=capselex$First_Optimal_offset[i], #12
                              optimal_orientation=capselex$First_Optimal_Orientation[i], #-
                              optimal_overlap=capselex$First_Optimal_Overlap[i], #8
                              dimer_length=metadata %>% filter(ID==capselex$first.monomer[i]) %>% pull(length)) #25
  
  
  second_monomer_se=monomer_start_end( tomtom_reverse_order=reverse_bool[[as.character(capselex$Second_tomtom_reverse_order[i])]], #FALSE
                     optimal_offset=capselex$Second_Optimal_offset[i], #-1
                     optimal_orientation=capselex$Second_Optimal_Orientation[i], #+#
                     optimal_overlap=capselex$Second_Optimal_Overlap[i], #9
                     dimer_length=metadata %>% filter(ID==capselex$second.monomer[i]) %>% pull(length))
  
  
  capselex$first_monomer_start[i]=first_monomer_se$start
  capselex$first_monomer_end[i]=first_monomer_se$end
  
  capselex$second_monomer_start[i]=second_monomer_se$start
  capselex$second_monomer_end[i]=second_monomer_se$end
  
  
  capselex$first_start[i]=first_se$start
  capselex$first_end[i]=first_se$end
  
  capselex$second_start[i]=second_se$start
  capselex$second_end[i]=second_se$end
  
  
  
  #cap length
  #Left and rightmost
  if(first_se$start <= second_se$start){
    
    gap=second_se$start-first_se$end #0
    
    if(gap>0){
      capselex$type_computation[i]="spacing"
      #if(gap==0){
      #  capselex$gap_length[i]=second_se$start-first_se$end
      #  capselex$gap_start[i]=first_se$end
      #  capselex$gap_end[i]=second_se$start
      #}else{
        capselex$gap_length[i]=second_se$start-first_se$end-1
        if(capselex$gap_length[i]!=0){
        capselex$gap_start[i]=first_se$end+1
        capselex$gap_end[i]=second_se$start
        }else{
          capselex$gap_start[i]=first_se$end
          capselex$gap_end[i]=first_se$end
        }
      #}
      
      
    }else{
      #gap is 0 or less
      capselex$type_computation[i]="composite"
      #The size of the overlapping part is abs(-3)+1=3
      capselex$overlap_length[i]=abs(gap)+1 #IS THIS CORRECT
      capselex$overlap_start[i]=second_se$start
      capselex$overlap_end[i]=first_se$end
      
    }
      
   }else{
     print("reverse")
     gap=first_se$start-second_se$end #0
     
     if(gap>0){
       capselex$type_computation[i]="spacing"
       #if(gap==0){
      #   capselex$gap_length[i]=first_se$start-second_se$end
      #   capselex$gap_start[i]=second_se$end
      #   capselex$gap_end[i]=first_se$start
         
      # }else{
         capselex$gap_length[i]=first_se$start-second_se$end-1
         if(capselex$gap_length[i]!=0){
          capselex$gap_start[i]=second_se$end+1
          capselex$gap_end[i]=first_se$start
         }else{
           capselex$gap_start[i]=second_se$end
           capselex$gap_end[i]=second_se$end
         }
       #}
     }else{
       capselex$type_computation[i]="composite"
       #The size of the overlapping part is abs(-3)+1=3
       capselex$overlap_length[i]=abs(gap)+1 #Is this correct
       capselex$overlap_start[i]=first_se$start
       capselex$overlap_end[i]=second_se$end
       
     }
  }
}


capselex$`type computation`=capselex$type_computation
capselex$overlap_length=-capselex$overlap_length

composite_ind=which(!is.na(capselex$overlap_length))
spacing_ind=which(!is.na(capselex$gap_length))

capselex$`type computation overlap or gap`=NA

capselex$`type computation overlap or gap`[composite_ind]=capselex$overlap_length[composite_ind]
capselex$`type computation overlap or gap`[spacing_ind]=capselex$gap_length[spacing_ind]

#How many of the new spacing motifs are predicted as composites
spacing_motifs=capselex %>% filter(study=="fromYimeng" & type=="spacing") #205

table(spacing_motifs$type_computation)
#          composite   spacing 
#12        37       156 

length(which(spacing_motifs$type_computation!="" & spacing_motifs$type==spacing_motifs$type_computation)) #156 CORRECT
length(which(spacing_motifs$type_computation!="" & spacing_motifs$type!=spacing_motifs$type_computation)) #37 UNCORRECT


#How many of the new composite motifs are predicted as spacing motifs
composite_motifs=capselex %>% filter(study=="fromYimeng" & type=="composite") #11131

length(which(composite_motifs$type_computation!="" & composite_motifs$type==composite_motifs$type_computation)) #1023 CORRECT
length(which(composite_motifs$type_computation!="" & composite_motifs$type!=composite_motifs$type_computation)) #93 UNCORRECT

table(composite_motifs$type_computation)
        #composite   spacing 
#15      1023        93 

#Reorder the first.monomer and second.monomer columns if switched_order==TRUE
#capselex=capselex[, c(1,2,25:50)]


#capselex=read_tsv(file="../../PWMs_final_union/metadata_google_drive_supplement_with_type_computation.tsv")

capselex$switched_order=FALSE
capselex$switched_order=capselex$first_start > capselex$second_start
table(capselex$switched_order)

#saveRDS(capselex,file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2.RDS")
saveRDS(capselex,file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2_artificial_halfsites.RDS")

# capselex_tomtom=capselex
# capselex_tomtom[which(capselex_tomtom$switched_order), c("first.monomer", "second.monomer")]=capselex_tomtom[which(capselex_tomtom$switched_order), c("second.monomer", "first.monomer")]
# capselex_tomtom=capselex_tomtom[,c("ID", "first.monomer", "second.monomer")]
# 
# matches=match(metadata$ID,capselex_tomtom$ID) 
# nona_ind=which(!is.na(matches))
# 
# capselex_tomtom$filename_dimer=metadata$filename[nona_ind]
# 
# matches=match(capselex_tomtom$first.monomer, metadata$ID) 
# nona_ind=which(!is.na(matches))
# 
# capselex_tomtom$filename_monomer1=metadata$filename[matches]
# 
# matches=match(capselex_tomtom$second.monomer, metadata$ID) 
# nona_ind=which(!is.na(matches))
# 
# capselex_tomtom$filename_monomer2=metadata$filename[matches]
# 
# write.table(capselex_tomtom, file="../../PWMs_final_version2.2/metadata_for_tomtom_figures.tsv",quote=FALSE, sep="\t", row.names = FALSE)







