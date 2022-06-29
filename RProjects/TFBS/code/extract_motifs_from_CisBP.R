

#############Read excel tables################

library(tidyverse)
library(readxl)


source("code/read_excel_tables.R")

# This table contains background subtracted PWMs for all of the TF pairs and individual TFs analyzed.  First line for each factor contains information as follows:	
# Base:	A, C, G or T
# symbol(s):	HGNC symbol (gene names)
# clone type:	Type of construct used: extended DNA binding domain (eDBD) or full-length (FL)
# families:	TF structural family
# experiment type:	Type of experiment: HT-SELEX or methyl-HT-SELEX
# ligand sequence:	Name of primer used to generate selection ligand; number before N indicates length of random sequence
# batch:	batch of HT-SELEX and Methyl-HT-SELEX experiment
# seed:	Sequence used to identify subsequences to be included in the model
# multinomial:	Hamming distance to seed used to identify subsequences to be included in the model
# cycle:	SELEX cycle used to construct the model; previous cycle was used as background unless indicated otherwise by a letter "b" followed by the used background cycle, "u" letter indicates that only the first instance of identical sequencing reads were used in the construction of the PWM
# Matrix is one of the representative PWMs:	value yes or no
# comment:	comment field


all_table <- read_excel("/home/osmalama/Dropbox/Taipale-lab/TFBS/PWMs/Nitta2015/CisBP_2022_06_22_4_11_am/PWMs.xlsx")

char_vector <- readLines("/home/osmalama/Dropbox/Taipale-lab/TFBS/PWMs/Nitta2015/CisBP_2022_06_22_4_11_am/PWM.txt")


tmp=which( diff( seq(1, length(char_vector), 1) [-which(char_vector=="")])==3)
tmp=as.list(tmp+seq(0,length(tmp)-1,1)*2)

for(i in 1:length(tmp)){
  if(i==1){
    start=1
  }else{
    start=last(tmp[[i-1]])+3
  }
  tmp[[i]]=seq(start,tmp[[i]],1)
  
  
}


# 
# 
# 

myfxn <- function(var1,var2){
  char_vector[var1]
}

PWMs_list=lapply(tmp, myfxn, var2=char_vector)


PWMs_metadata=data.frame(TF=character(),
           TF_Name=character(),
           Gene=character(),
           Motif=character(),
           Family=character(),
           Species=character(),
           Motif_missing=character())

dir.create(paste0("../../PWMs/Nitta2015/CisBP_2022_06_22_4_11_am/pwms/","Homo_sapiens"), recursive=TRUE)

for(m in 1:length(PWMs_list)){
  PWMs_metadata[m,]=c( do.call(rbind, strsplit(PWMs_list[[m]][1:6],split="\t"))[,2], "No")
 
  #The motif
  Motif=do.call(rbind, strsplit(PWMs_list[[m]][7:length(PWMs_list[[m]])],split="\t"))
  if(nrow(Motif)==1){
    PWMs_metadata[m,7]="Yes"
  }else{
    #Write motif
    #A C G T
    write.table(t(Motif[-1, -1]),file=paste0("../../PWMs/Nitta2015/CisBP_2022_06_22_4_11_am/pwms/","Homo_sapiens/", 
                                                         paste0(PWMs_metadata[m,
                                                                              -which(colnames(PWMs_metadata)%in% c("Gene", "Family", "Species", "Motif_missing"))], collapse="_"),
                                                         ".pwm"),
                quote=FALSE, row.names=FALSE, col.names=FALSE
                                                         
                                                         )
    
  }
   
}




colnames(PWMs_metadata)=c("symbol",	"clone type","families",	"experiment",
                          "ligand sequence",	"batch", "seed",	"multinomial",
                          "cycle","Matrix is one of the representative PWMs", "comment")

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone type","families",	"experiment",
                                     "ligand sequence",	"batch", "seed",	"multinomial",
                                     "cycle","Matrix is one of the representative PWMs", "comment")
)
PWMs_metadata$organism="Homo_sapiens"

PWMs_metadata=PWMs_metadata %>% 
  mutate(families = gsub("/|\\.\\d+[A-Za-z]+", "_aka_", families))



PWMs_metadata %>%
  count(symbol) 


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

dir.create(paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens"), recursive=TRUE)
for(m in 1:length(PWMs_list)){
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0("../../PWMs/Yin2017/pwms/","Homo_sapiens","/", 
                          paste0(PWMs_metadata[m,
                          -which(colnames(PWMs_metadata)%in% c("organism", "comment"))], collapse="_"),
                          ".pfm"))
}





#display_table_shape(all_table)

