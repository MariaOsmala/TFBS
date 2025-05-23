

#############Read excel tables################

library(tidyverse)
library(readxl)
#all_table <- read_excel("useful/tidyverse_notes/utility/multiple_tables_sheet.xlsx", col_names=FALSE)

#  his table contains background subtracted PWMs for all factors analyzed.  First line for each factor contains information as follows:	
# symbol:	HGNC symbol (gene name); Human is all caps, in mouse constructs, only the first initial is capitalized
# family:	TF structural family
# clone type:	Type of construct used: DBD or full-length
# ligand sequence:	Name of primer used to generate selection ligand; number before N indicates length of random sequence
# batch:	batch of SELEX experiment
# N-length:	Length of randomized region
# seed:	Sequence used to identify subsequences to be included in the model
# multinomial:	Hamming distance to seed used to identify subsequences to be included in the model
# cycle:	SELEX cycle used to construct the model; previous cycle was used as background
# site: type	number of repeats of a subsequence (monomer, dimer, tetramer or dimer of dimers)
# comment:	comment field
# Matrix is  one of the representative PWMs	value: yes or no

# Jolma2013
# Models indicated in orange were removed from network representation of data (Figure 3 and supplementary figures) 
# due to reason indicated in the comment field. Models marked green are technical replicates used in counting of error bar in figure 1D.	
# 
oranges=which(seq(18, 4228,5) %in% seq(4118,4168,5)) #11
greens=which(seq(18, 4228,5) %in% seq(4173,4228,5)) #12

all_table <- read_excel("/home/osmalama/Dropbox/Taipale-lab/TFBS/PWMs/Jolma2013/mmc3.xls", col_names=FALSE, skip=17)

PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,]))

#Jolma2013
colnames(PWMs_metadata)=c("symbol",	"family",	"clone type",
                          "ligand sequence",	"batch", "Seed",	"multinomial",
                          "cycle", "site type",	"comment","Matrix is  one of the representative PWMs")


PWMs_metadata=select(PWMs_metadata, c("symbol",	"family",	"clone type",
                                      "ligand sequence",	"batch", "Seed",	"multinomial",
                                      "cycle", "site type",	"comment","Matrix is  one of the representative PWMs"))



#Jolma2015



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



all_table <- read_excel("/home/osmalama/Dropbox/Taipale-lab/TFBS/PWMs/Jolma2015/41586_2015_BFnature15518_MOESM33_ESM.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=19)

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,-1]))

colnames(PWMs_metadata)=c("symbol",	"family",	"experiment",
                          "ligand sequence",	"batch", "seed",	"multinomial",
                          "cycle","Matrix is one of the representative PWMs", "comment",
                          "conservation p-value", "ChIP-seq cluster p-value", "category")

PWMs_metadata=select(PWMs_metadata,c("symbol",	"family",	"experiment",
                                     "ligand sequence",	"batch", "seed",	"multinomial",
                                     "cycle","Matrix is one of the representative PWMs","comment", 
                                     "conservation p-value", "ChIP-seq cluster p-value", "category"))
PWMs_metadata$organism="Homo_sapiens"


PWMs_metadata %>%
  count(symbol) 
#463 unique

##Add organism info

PWMs_metadata =PWMs_metadata %>%
  mutate(
   organism = if_else(
      condition = str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]'), 
      true      = "Homo_sapiens", 
      false     = "Mus_musculus"
    )
  )

#Number of human TFs
PWMs_metadata %>%
  count(symbol)%>%
  filter(str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]'))
#381   

#Number of human TF families
PWMs_metadata %>%
  filter(str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]')) %>%
    count(family)

#38   


#Number of mouse TFs
PWMs_metadata %>%
  count(symbol)%>%
  filter(str_detect(symbol, '[:lower:]'))
#82

#Number of mouse TF families
PWMs_metadata %>%
  filter(str_detect(symbol, '[:lower:]')) %>%
  count(family)
#13   


#filter oranges and greens
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  count(symbol)
#463

#Number of Human TFs
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  filter(str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]'))%>%
  count(symbol)
#381

#Number of Human TFfamilies
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  filter(str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]'))%>%
  count(family)
#37

#Number of mouse TFs
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  count(symbol)%>%
  filter(str_detect(symbol, '[:lower:]'))
#82

#Number of mouse TFs
PWMs_metadata %>%
  slice(-c(oranges, greens))%>%
  filter(str_detect(symbol, '[:lower:]'))%>%
  count(family)
#13

#How many types unique(PWMs_metadata$`site type`)
#[1] "monomeric"             "dimeric"               "monomer"               "monomeric or dimeric"  "dimer of dimers"       "putative multimer"    
#[7] "monomer or dimer"      "trimeric"              "putatively multimeric"

#Number of human TFs representative
PWMs_metadata %>%
  filter(str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]')) %>%
  filter(`Matrix is  one of the representative PWMs`=="no")%>%
  count(symbol)

#representative 137
#nonrepresentative 303

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Jolma2013/pwms/",PWMs_metadata[m,"organism"],"/", paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("comment", "organism"))], collapse="_"),".pfm"))
}





#display_table_shape(all_table)

