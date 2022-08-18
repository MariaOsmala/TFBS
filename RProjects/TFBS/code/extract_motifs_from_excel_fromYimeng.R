

#############Read excel tables################

library(tidyverse)

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

types=c("pfm_composite_new", "pfm_spacing_new")

path="../../PWMs/fromYimeng/pfm_from_newData/"

dir("../../PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/")

datas=tibble()
  
for(t in types){
  data=as.tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),"_")))
  data$V8=gsub(".pfm", "",as.character(data$V8))
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$study="fromYimeng"
  data$type=t
  data$organism="Homo_sapiens"
  data$clone=NA
  data$family=NA
  data$experiment="CAP-SELEX"
  data$representative=NA
  data$comment=NA
  data$filename=paste0("PWMs/fromYimeng/pfm_from_newData/", t,"/", dir(paste0(path, t, "/")))
  datas=rbind(datas, data)
}

datas=datas[, c("symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","representative", "short","type","comment","filename")]

datas$multinomial=gsub("m", "",as.character(datas$multinomial))
datas$cycle=gsub("c", "",as.character(datas$cycle))

write.table(datas, file="../../PWMs/fromYimeng/metadata.csv", row.names = FALSE)
saveRDS(datas, file="data/fromYimeng.Rds")


#display_table_shape(all_table)

