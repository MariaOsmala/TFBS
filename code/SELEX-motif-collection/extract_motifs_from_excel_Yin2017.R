

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
#library(DECIPHER)
library(Biostrings)
library(TFBSTools)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)


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


all_table <- read_excel("../../Data/SELEX-motif-collection/aaj2239_yin_sm_tables_s1-s6.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=22, n_max=8996-22, col_types=c(rep("text",12), rep("numeric",19)))

 
PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(seq(1,nrow(PWMs),5), function(i) PWMs[i,-1]))



colnames(PWMs_metadata)=c("symbol",	"clone","family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment")
PWMs_metadata$organism="Homo_sapiens"
PWMs_metadata$study="Yin2017"
PWMs_metadata$short=NA
PWMs_metadata$type=NA


#symbol     clone family organism     study      experiment ligand        batch  seed               multinomial cycle representative short type              comment filename         

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
                                     "ligand",	"batch", "seed",	"multinomial",
                                     "cycle","representative", "short", "type","comment")
)


PWMs_metadata=PWMs_metadata %>% 
  mutate(family = gsub("/|\\.\\d+[A-Za-z]+", "_aka_", family))


PWMs_metadata$filename=NA
PWMs_metadata$ID=NA
PWMs_metadata$IC=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=NA


PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) ) #1794


append=FALSE



#Remove these from metadata
remove=c()

for(m in 1:length(PWMs_list)){
  
  #Empty matrix
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)    
  #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    
  }else{
  
    
  # Write .pfm tab-separated
  filename=paste0(pfms_tab_path, "/",
                  paste0(PWMs_metadata[m,
                                       -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                            "representative", "short", "type", "filename","ID", "IC", "length", "consensus"))], collapse="_"),
                  ".pfm")
  
  ID=paste0(PWMs_metadata[m,
                       -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                            "ID", "IC", "length", "consensus" ))], collapse="_")
  
  
  #Filename and ID should be unique
  if(filename %in% PWMs_metadata$filename){
    print("Warning! Non-unique filenames and IDs")
  }
  
  
  PWMs_metadata$filename[m]=filename
  PWMs_metadata$ID[m]=ID
  
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=filename, sep="\t")
  
  pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
  colnames(pcm)=NULL
  rownames(pcm)=c("A", "C", "G", "T")
  
  pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                              profileMatrix=pcm
  )
  
  
  pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
  PWMs_metadata[m, "IC"]=sum(rowSums(icm)) 
  
  #motif length
  PWMs_metadata[m, "length"]=length(pfm)
  
  motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
  motif@name=paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "representative","short", "type", "filename",
                                                                         "ID", "IC", "length", "consensus"))], collapse="_")
  
  PWMs_metadata[m, "consensus"]=motif@consensus
  
  #Write transfac format
  transfac=paste0(pfms_transfac_path,"/",
                  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone","family","comment", "study","organism","representative","short", "type","filename",
                                                                              "ID", "IC", "length", "consensus"))], collapse="_"),".pfm")
  write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
  #Write .pfm space separated
  write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_space_path, "/",
                          paste0(PWMs_metadata[m,
                                               -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative","short", "type", "filename",
                                                                                    "ID", "IC", "length", "consensus"))], collapse="_"),
                          ".pfm"), sep=" ")
  
  #Write .scpd format, all into a single file
  PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
  rownames(PWM)=c("A", "C", "G", "T")
  write.table(paste0(">",  paste0(PWMs_metadata[m,-which(colnames(PWMs_metadata)%in% c("clone", "family", "organism", "study","comment", "representative", "short", "type", "filename",
                                                                                       "ID", "IC", "length", "consensus"))], collapse="_")),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/", "Yin2017_all", ".scpd"))
  append=TRUE
  
  write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/", "Yin2017_all", ".scpd"))
  
  
  #Write also pwms 
  
  write.table(pwm_class@profileMatrix, file=paste0(pwms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
  write.table(pwm_class@profileMatrix, file=paste0(pwms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
  
  
  # Draw logos --------------------------------------------------------------
  
  
  
  png(paste0(logos_png_prob,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
          ggseqlogo::theme_logo()  +
          labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  
  dev.off()
  
  
  pdf(paste0(logos_pdf_prob,"/",ID,".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
          ggseqlogo::theme_logo()  +
          labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  dev.off()
  
  #Information content matrices
  
  
  png(paste0(logos_png_ic,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
          ggseqlogo::theme_logo()  +
          labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  
  dev.off()
  
  
  pdf(paste0(logos_pdf_ic,"/","ID",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
          ggseqlogo::theme_logo()  +
          labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  dev.off()
  
  
  
  }
}

PWMs_metadata=PWMs_metadata[-remove,]

#Add Family Info



#HumanTFs-ccbr-TableS1.csv
TF_family_mappings <- read_csv("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol","Gene Information/ID","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)
PWMs_metadata$DBD=NULL


PWMs_metadata <- PWMs_metadata %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"= "Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")



Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") #923

Methyl_SELEX_motif_categories=read_delim("../../Data/SELEX-motif-collection/Methyl_SELEX_metadata_with_info.csv", 
                                                                            delim = ";", 
                                         escape_double = FALSE, trim_ws = TRUE)

PWMs_metadata=PWMs_metadata %>%
  left_join(select(Methyl_SELEX_motif_categories, "ID", "final suggestion" ), by = 'ID')

colnames(PWMs_metadata)[23]="Methyl.SELEX.Motif.Category"

table(PWMs_metadata$Methyl.SELEX.Motif.Category) #some typos

PWMs_metadata$Methyl.SELEX.Motif.Category[PWMs_metadata$Methyl.SELEX.Motif.Category=="No GpG"]="No CpG"
PWMs_metadata$Methyl.SELEX.Motif.Category[PWMs_metadata$Methyl.SELEX.Motif.Category=="Multiple Effects"]="Multiple effects"

categories=table(PWMs_metadata$Methyl.SELEX.Motif.Category)

pie(categories, 
    main = "Proportions of different Methyl-SELEX motif categories", 
    col = rainbow(length(categories)), 
    labels = paste(names(categories), "\n", categories))

remove=c("Little effect", "MethylMinus", "No CpG")

PWMs_metadata=PWMs_metadata[-which(PWMs_metadata$Methyl.SELEX.Motif.Category %in% remove),]

sum(table(PWMs_metadata$Methyl.SELEX.Motif.Category)) #297

#HT-SELEX Methyl-HT-SELEX 
#864             297 

write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Yin2017_metadata.csv", row.names = FALSE,sep="\t")


