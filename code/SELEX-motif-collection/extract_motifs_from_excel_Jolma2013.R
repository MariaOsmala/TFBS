

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
#all_table <- read_excel("useful/tidyverse_notes/utility/multiple_tables_sheet.xlsx", col_names=FALSE)

source("code/read_excel_tables.R")
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
#Remove these from the models
oranges=which(seq(18, 4228,5) %in% seq(4118,4168,5)) #11
greens=which(seq(18, 4228,5) %in% seq(4173,4228,5)) #12

all_table <- read_excel("../../PWMs/Jolma2013/mmc3.xls", sheet="Table S3 - PWM models", col_names=FALSE, skip=17, col_types=c(rep("text",11), rep("numeric",13)))

PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,]))

#Jolma2013
colnames(PWMs_metadata)=c("symbol",	"family",	"clone",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle", "type",	"comment","representative")

PWMs_metadata=select(PWMs_metadata, c("symbol",	"family",	"clone",
                                      "ligand",	"batch", "seed",	"multinomial",
                                      "cycle", "type",	"comment","representative"))


##Add organism info

PWMs_metadata =PWMs_metadata %>%
  mutate(
    organism = if_else(
      condition = str_detect(symbol, '[:upper:]') & !str_detect(symbol, '[:lower:]'), 
      true      = "Homo_sapiens", 
      false     = "Mus_musculus"
    )
  )

PWMs_metadata$study="Jolma2013"
PWMs_metadata$experiment="HT-SELEX"
PWMs_metadata$short=NA
PWMs_metadata$filename=NA


PWMs_metadata$representative=toupper(PWMs_metadata$representative)

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment", "filename"))


PWMs_metadata$ID=""
PWMs_metadata$IC=NA
PWMs_metadata$IC_universal=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=""

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

PWMs_metadata=PWMs_metadata[-c(oranges, greens),]
PWMs_list=PWMs_list[-c(oranges, greens)]

pwms_human="../../PWMs_final/Jolma2013/pwms/Homo_sapiens"
pwms_mouse="../../PWMs_final/Jolma2013/pwms/Mus_musculus"

pwms_space_human="../../PWMs_final/Jolma2013/pwms_space/Homo_sapiens"
pwms_space_mouse="../../PWMs_final/Jolma2013/pwms_space/Mus_musculus"

transfac_human="../../PWMs_final/Jolma2013/transfac/Homo_sapiens"
transfac_mouse="../../PWMs_final/Jolma2013/transfac/Mus_musculus"

dir.create(pwms_human, recursive=TRUE)
dir.create(pwms_space_human, recursive=TRUE)
dir.create(transfac_human, recursive=TRUE)
dir.create(pwms_mouse, recursive=TRUE)
dir.create(pwms_space_mouse, recursive=TRUE)
dir.create(transfac_mouse, recursive=TRUE)



#Remove these from metadata
remove=c()
pwms=list("Homo_sapiens"=pwms_human, "Mus_musculus"=pwms_mouse)
pwms_space=list("Homo_sapiens"=pwms_space_human, "Mus_musculus"=pwms_space_mouse)
transfac_path=list("Homo_sapiens"=transfac_human, "Mus_musculus"=transfac_mouse)

append=FALSE 



for(m in 1:length(PWMs_list)){
  
  
  #Empty matrix
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
  }else{
    
  # Write .pfm tab-separated
  filename=paste0(pwms[[PWMs_metadata$organism[m]]],"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                              "representative", "short", "type", "filename","ID", "IC", "IC_universal","length", "consensus"))], collapse="_"),
                    ".pfm")
  ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                 "ID", "IC", "IC_universal","length", "consensus" ))], collapse="_")
  #Filename and ID should be unique
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
    
  
    #Write pfm
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE, file=filename,sep="\t")
  
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
     
    transfac=paste0(transfac_path[[ PWMs_metadata$organism[m] ]],"/",ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  

    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0( pwms_space[[PWMs_metadata$organism[m]]],"/",ID,
                            ".pfm"),sep=" ")
  
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Jolma2013/all", ".scpd"))
    append=TRUE
  
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs_final/Jolma2013/","all", ".scpd"))
  
  }
}

#convert to .scpd format (this can be converted to meme format)
# >BCL6B_1
# A	0	3	0	0	0	0	0	0	280	0	0	289	290	0	0	113	183
# C	19	0	367	0	0	0	384	12	0	0	0	0	0	2	16	243	146
# G	0	186	0	0	0	2	0	34	0	173	171	2	0	0	0	9	7
# T	207	1	2	212	213	212	0	214	0	0	0	0	0	215	211	0	0

TF_family_mappings <- read_csv("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)

PWMs_metadata$`Gene Information/ID`=NULL
PWMs_metadata$DBD=NULL

write.table(PWMs_metadata, file="../../PWMs_final/Jolma2013/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(PWMs_metadata, file="Rdata/Jolma2013.Rds")


stop()




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






#display_table_shape(all_table)

