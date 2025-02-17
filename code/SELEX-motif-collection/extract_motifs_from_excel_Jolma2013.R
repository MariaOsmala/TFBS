

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
config <- yaml.load_file("../../Experiments/config.yaml", readLines.warn=FALSE)

source("paths_to_folders.R", echo=FALSE)

# This table contains background subtracted PWMs for all factors analyzed.  First line for each factor contains information as follows:	
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
# Matrix is one of the representative PWMs	value: yes or no

# Jolma2013
# Models indicated in orange were removed from network representation of data (Figure 3 and supplementary figures) 
# due to reason indicated in the comment field. Models marked green are technical replicates used in counting of error bar in figure 1D.	
# 
#Remove these from the models
oranges=which(seq(18, 4228,5) %in% seq(4118,4168,5)) #11
greens=which(seq(18, 4228,5) %in% seq(4173,4228,5)) #12

all_table <- read_excel("../../Data/SELEX-motif-collection/mmc3.xls", sheet="Table S3 - PWM models", col_names=FALSE, skip=17, col_types=c(rep("text",11), rep("numeric",13)))
#display_table_shape(all_table)

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
PWMs_metadata$length=NA
PWMs_metadata$consensus=""

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

PWMs_metadata=PWMs_metadata[-c(oranges, greens),]
PWMs_list=PWMs_list[-c(oranges, greens)]




#Remove these from metadata
remove=c()

append=FALSE 



for(m in 1:length(PWMs_list)){ #820
  
  #Remove if just empty? matrix
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
  }else{
    
  # Write .pfm tab-separated
  filename=paste0(pfms_tab_path,"/", 
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
    
  
    #Write pfm
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE, file=filename,sep="\t")
  
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
    
    
    #Write transfac
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    transfac=paste0(pfms_transfac_path,"/",ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
    #consensus obtained from universalmotif
    PWMs_metadata[m, "consensus"]=motif@consensus
    
    
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0( pfms_space_path,"/",ID,
                            ".pfm"),sep=" ")
  
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","Jolma2013", ".scpd"))
    append=TRUE
  
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","Jolma2013", ".scpd"))
    
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



TF_family_mappings <- read_csv("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv")


# Merge two data frames
PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID","DBD" ), by = 'symbol')

PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = PWMs_metadata$DBD, .before = 4)

#Rename Gene Information/ID column and relocate column 
#PWMs_metadata <- PWMs_metadata %>%  
#  rename("new_name" = "old_name")


PWMs_metadata <- PWMs_metadata %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"= "Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

PWMs_metadata$DBD=NULL

write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Jolma2013_metadata.csv", row.names = FALSE, sep="\t")








