library("universalmotif") #does not work

library("devtools")
install_github("Danko-Lab/rtfbs_db/rtfbsdb")

library(rtfbsdb);
#Create db object from pre-installed dataset
db <- CisBP.extdata("Homo_sapiens");
db <- CisBP.extdata("Mus_musculus");
db <- CisBP.extdata("Drosophila_melanogaster");

# Create db object from a zipped downloaded dataset
db <- CisBP.zipload("../../PWMs/Jolma2013/CisBP_2022_06_22_4_02_am/data.zip", "Homo_sapiens");
tfs <- tfbs.createFromCisBP(db);
# Create db object from downloaded dataset
db <- CisBP.download("Mus_musculus");



cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",package = "universalmotif"))
write_matrix(cisbp, file="../../PWMs/tmp.txt")

#Read motifs, convert to MOODS format
PWM.list=read_cisbp("../../PWMs/Jolma2013/CisBP_2022_06_22_4_02_am/data.zip", skip=0)


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
# Models indicated in orange were removed from network representation of data (Figure 3 and supplementary figures) 
# due to reason indicated in the comment field. Models marked green are technical replicates used in counting of error bar in figure 1D.	
# 

oranges=which(seq(18, 4228,5) %in% seq(4118,4168,5)) #11

greens=which(seq(18, 4228,5) %in% seq(4173,4228,5)) #12



green

all_table <- read_excel("/home/osmalama/Dropbox/Taipale-lab/TFBS/PWMs/Jolma2013/mmc3.xls", col_names=FALSE, skip=17)
  
 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,]))

colnames(PWMs_metadata)=c("symbol",	"family",	"clone type",
                          "ligand sequence",	"batch", "Seed",	"multinomial",
                          "cycle", "site type",	"comment","Matrix is  one of the representative PWMs")


PWMs_metadata=select(PWMs_metadata, c("symbol",	"family",	"clone type",
                                      "ligand sequence",	"batch", "Seed",	"multinomial",
                                      "cycle", "site type",	"comment","Matrix is  one of the representative PWMs"))


# PWM.summary
# Exp                           Species Unique.n TFs.n Family.n
# 2 Jolma2013 Homo_sapiens      372   372       36
# 3 Jolma2013 Mus_musculus       80    80       12
# colSums(PWM.summary[2:3, 3:5])
# Unique.n    TFs.n Family.n 
# 452      452       48 
 

PWMs_metadata %>%
  count(symbol) 
#463 unique

#Add organism info

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
  filter(`Matrix is  one of the representative PWMs`=="yes")%>%
  count(symbol)



PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Jolma2013/pwms/",PWMs_metadata[m,"organism"],"/", paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("comment", "organism"))], collapse="_"),".pfm"))
}

#metadata in the first rows




#display_table_shape(all_table)

