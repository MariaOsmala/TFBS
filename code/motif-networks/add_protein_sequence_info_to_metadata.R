library("readr")
library("tidyverse")

rm(list=ls())



load("~/projects/SELEX/TF-amino-acid-sequences/RData/protein_sequences.RData")

metadata=read.table("~/projects/TFBS/code/motif-networks/motif_webpage/metadata.tsv",
                    sep="\t", header=TRUE, stringsAsFactors = FALSE)

metadata$logo_path=NULL


metadata$`protein sequence 1`=NA
metadata$`protein sequence 2`=NA


#What info we have for protein sequences 
setdiff(names(df_Jolma2013_protein), names(metadata))
# "HNGC" "Ensembl_ID" "protein_sequence" "DOMAIN" "source organism" 
# "clone type" "production method"                       
# "Clone source" "Classificaton in Vaquerizas et al., 2009"



setdiff(names(df_Jolma2015_monomers_protein), names(metadata))
#"HNGC-name"                   "Ensembl id."                 "amino-acid sequence"        
#"main structural class"       "6-mer"                       "source organism"            
#"Clone source"                "Construct type"             "Included domain Ids"        
#"Included domain description" "Excluded domain Ids"         "Excluded domain description"

#Remove these columns
#"Protein amount"             
#"HT-SELEX validation"         "Activity in CAP-SELEX"      


setdiff(names(df_Jolma2015_heterodimers_protein), names(metadata))

# "symbol_both"                     "TF1"                            
# "TF1 HNGC-name"                   "TF1 Ensembl id."                
# "TF1 amino-acid sequence"         "TF1 main structural class"      
# "TF1 6-mer"                       "TF1 source organism"            
# "TF1 Clone source"                "TF1 Construct type"             
# "TF1 Protein amount"              "TF1 HT-SELEX validation"        
# "TF1 Activity in CAP-SELEX"       "TF1 Included domain Ids"        
# "TF1 Included domain description" "TF1 Excluded domain Ids"        
# "TF1 Excluded domain description" "TF2 HNGC-name"                  
# "TF2 Ensembl id."                 "TF2 amino-acid sequence"        
# "TF2 main structural class"       "TF2 6-mer"                      
# "TF2 source organism"             "TF2 Clone source"               
# "TF2 Construct type"              "TF2 Protein amount"             
# "TF2 HT-SELEX validation"         "TF2 Activity in CAP-SELEX"      
# "TF2 Included domain Ids"         "TF2 Included domain description"
# "TF2 Excluded domain Ids"         "TF2 Excluded domain description"


setdiff(names(df_Yin2017_protein), names(metadata))
#"symbol_clone"                             "Ensembl_ID"                               "HNGC"                                    
#"clone type"                               "DOMAIN"                                   "Classificaton in Vaquerizas et al., 2009"
#"protein_sequence"                         "source organism"                          "production method"                       
#"Clone source"     

setdiff(names(df_Xie2025_monomers_protein), names(metadata))
#"Ensembl_ID"    "HNGC"          "Clone used"    "DOMAIN"        "eDBD seuqence"

setdiff(names(df_Xie2025_heterodimers_protein), names(metadata))
#"symbol_both"       "TF1"               "TF1 Ensembl_ID."   "TF1 HNGC"          "TF1 Clone used"    "TF1 DOMAIN"        "TF1 eDBD sequence"
#"TF2 Ensembl_ID."   "TF2 HNGC"          "TF2 Clone used"    "TF2 DOMAIN"        "TF2 eDBD sequence"


# Jolma2013 ---------------------------------------------------------------

df_Jolma2013_protein$type %>% table(useNA="always")
#dimer of dimers               dimeric               monomer      monomer or dimer             monomeric 
#7                             434                   1            1                            355 
#monomeric or dimeric     putative multimer putatively multimeric              trimeric                  <NA> 
#1                        5                 9                                  7                     0 

df_Jolma2013_protein %>% filter(type %in% c("monomer", "monomeric", "monomer or dimer", "monomeric or dimeric")) 


metadata=metadata %>% #820 x 25
  left_join(df_Jolma2013_protein %>% 
              select(ID, "protein_sequence"),
            by="ID") %>%
  mutate(`protein sequence 1` = protein_sequence) %>%
  select(-protein_sequence)


# Yin2017 -----------------------------------------------------------------

names(table(df_Yin2017_protein$ID[grep("ensembl ID is different", df_Yin2017_protein$`Classificaton in Vaquerizas et al., 2009`)]))

grep("ensembl ID is different", df_Yin2017_protein$`Classificaton in Vaquerizas et al., 2009`) 
#102 103 226 378 404 447 543 544 802 832

which(df_Yin2017_protein$ID %in% names((df_Yin2017_protein$ID %>% table())[which((df_Yin2017_protein$ID %>% table())>1 )]))
#377 378 542 543 801 802

#remove 378, 543, 802, 

df_Yin2017_protein=df_Yin2017_protein[-c(378,543,802),]


Yin2017_data<- df_Yin2017_protein %>%
  select(ID, protein_sequence)

#Do not know what these are so add only once

metadata<- metadata %>%
  left_join(Yin2017_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`)
  ) %>%
  select(-protein_sequence)


Yin2017_Methyl_data<- df_Yin2017_Methyl_protein %>%
  select(ID, protein_sequence)

#Do not know what these are so add only once

metadata<- metadata %>%
  left_join(Yin2017_Methyl_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(protein_sequence), protein_sequence, `protein sequence 1`)
  ) %>%
  select(-protein_sequence)



# Jolma2015 ---------------------------------------------------------------

#Did not find SELEX data for monomers
Jolma2015_data<- df_Jolma2015_monomers_protein %>%
  select(ID, `amino-acid sequence`)

metadata <- metadata %>%
  left_join(Jolma2015_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`amino-acid sequence`), `amino-acid sequence`, `protein sequence 1`)
  ) %>%
  select(-`amino-acid sequence`)


Jolma2015_data<- df_Jolma2015_heterodimers_protein %>%
  select(ID, `TF1 amino-acid sequence`, `TF2 amino-acid sequence`)


metadata <- metadata %>%
  left_join(Jolma2015_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`TF1 amino-acid sequence`), `TF1 amino-acid sequence`, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(`TF2 amino-acid sequence`), `TF2 amino-acid sequence`, `protein sequence 2`)
  ) %>%
  select(-c(`TF1 amino-acid sequence`, `TF2 amino-acid sequence`))


# Xie2025 -----------------------------------------------------------------


#Did not find SELEX data for monomers
Xie2025_data<- df_Xie2025_monomers_protein %>%
  select(ID, `eDBD seuqence`)

metadata <- metadata %>%
  left_join(Xie2025_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`eDBD seuqence`), `eDBD seuqence`, `protein sequence 1`)
  ) %>%
  select(-`eDBD seuqence`)


Xie2025_data<- df_Xie2025_heterodimers_protein %>%
  select(ID, `TF1 eDBD sequence`, `TF2 eDBD sequence`)


metadata <- metadata%>%
  left_join(Xie2025_data, by = "ID") %>%
  mutate(
    `protein sequence 1` = if_else(!is.na(`TF1 eDBD sequence`), `TF1 eDBD sequence`, `protein sequence 1`),
    `protein sequence 2` = if_else(!is.na(`TF2 eDBD sequence`), `TF2 eDBD sequence`, `protein sequence 2`)
  ) %>%
  select(-c(`TF1 eDBD sequence`, `TF2 eDBD sequence`))


# Stats -------------------------------------------------------------------

#For how many we were able to find protein sequences

nrow(metadata) #3933

is.na(metadata$`protein sequence 1`) %>% table()

metadata%>% filter( is.na(`protein sequence 1`)) %>% pull(ID)

write_delim(metadata, "~/projects/TFBS/code/motif-networks/motif_webpage/metadata_with_protein_sequences.tsv", delim="\t")

