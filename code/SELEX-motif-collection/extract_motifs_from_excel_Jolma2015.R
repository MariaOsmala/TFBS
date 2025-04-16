

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)
source("compute_kl_divergence.R", echo=FALSE)
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

oranges=which(seq(20, 3335,5) %in% seq(2830,3160,5)) #67
lilas=which(seq(20, 3335,5) %in% seq(3165,3325,5)) #33
browns=which(seq(20, 3335,5) %in% seq(3330,3335,5)) #2

all_table <- read_excel("../../Data/SELEX-motif-collection/41586_2015_BFnature15518_MOESM33_ESM.xlsx", 
                        sheet="S2 PWM models", col_names=FALSE, skip=19, col_types=c(rep("text",13), rep("numeric",20)))

PWMs=split_df(all_table)[[1]]

PWMs_metadata=do.call(rbind,lapply(
                                  seq(1,nrow(PWMs),5), 
                                  function(i) PWMs[i,-1]
                                  )
                      )







colnames(PWMs_metadata)=c("symbol",	"family",	"experiment",
                          "ligand",	"batch", "seed",	"multinomial",
                          "cycle","representative", "comment",
                          "conservation p-value", "ChIP-seq cluster p-value", "category")

#Remove Hoxa10 Pantheropsis Guttatus, snake
snake_ind=which(PWMs_metadata$symbol=="Hoxa10")


#remove ChIP-EXO motifs
PWMs_metadata=PWMs_metadata[-c(snake_ind,oranges,browns),]

#"symbol",	"clone","family", "organism",	"study","experiment",
#"ligand",	"batch", "seed",	"multinomial",
#"cycle","representative", "short", "type","comment", "filename"

PWMs_metadata$clone=NA
PWMs_metadata$study="Jolma2015"
PWMs_metadata$short=NA
PWMs_metadata$type=NA
PWMs_metadata$filename=NA
PWMs_metadata$organism="Homo_sapiens"

PWMs_metadata$representative=toupper(PWMs_metadata$representative)

PWMs_metadata=select(PWMs_metadata,c("symbol",	"clone","family", "organism",	"study","experiment",
  "ligand",	"batch", "seed",	"multinomial",
  "cycle","representative", "short", "type","comment", "filename"))


PWMs_metadata$ID=""
PWMs_metadata$IC=NA
PWMs_metadata$length=NA
PWMs_metadata$consensus=""
PWMs_metadata$kld_between_revcomp=NA #Kullback-Leibler divergence between the motif and its reverse complement. Gives idea whether the motif is palindromic.

#no representative info
ind=which(!PWMs_metadata$representative %in% c("YES", "NO")) # lilas and ind are the same
PWMs_metadata$representative[ind]=NA

PWMs_list=lapply(seq(1,nrow(PWMs),5), function(i) PWMs[(i+1):(i+4),] %>%
                   select_if(~ !any(is.na(.))) )

#remove ChIP-EXO motifs and replicates
PWMs_list=PWMs_list[-c(snake_ind,oranges,browns)]

append=FALSE


remove=c()



for(m in 1:length(PWMs_list)){
  
  numeric_matrix <- matrix(as.numeric(as.vector(as.matrix(PWMs_list[[m ]][,-1]))), nrow = nrow(as.matrix(PWMs_list[[m ]][,-1])))
  
  if( ncol(PWMs_list[[m ]]) ==1 ) {
    remove=c(remove, m)
    print(m)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWMs_list[[ m ]][,-1]))) ==0 ){
    remove=c(remove, m)    
    print(m)
   
   }else{
     #There are zero columns in the matrix, remove them
    if(length(which(colSums(numeric_matrix)==0))>0){
      # Assume 'char_matrix' is your character matrix
      numeric_matrix <- matrix(as.numeric(as.vector(as.matrix(PWMs_list[[m ]][,-1]))), nrow = nrow(as.matrix(PWMs_list[[m ]][,-1])))
      print(m)
      #print(PWMs_list[[m ]])
      #Remove the zero column, it is always the last
      PWMs_list[[m ]]=PWMs_list[[m ]][,-ncol(PWMs_list[[m ]])]
      
    }
    
    
    
    
    # Write .pfm tab-separated
    filename=paste0(pfms_tab_path,"/", 
                    paste0(PWMs_metadata[m,
                                         -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                              "representative", "short", "type", "filename","ID", "IC","length", "consensus","kld_between_revcomp"))], collapse="_"),
                    ".pfm")
    ID=paste0(PWMs_metadata[m,
                            -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                 "ID", "IC", "length", "consensus","kld_between_revcomp" ))], collapse="_")
    
    if(filename %in% PWMs_metadata$filename){
      print("Warning! Non-unique filenames and IDs")
      print(ID)
      filename=paste0(pfms_tab_path,"/", 
                      paste0(PWMs_metadata[m,
                                           -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", 
                                                                                "representative", "short", "type", "filename","ID", "IC", "length", "consensus","kld_between_revcomp"))], collapse="_"),
                      "_v2",".pfm")
      ID=paste0( paste0(PWMs_metadata[m,
                              -which(colnames(PWMs_metadata)%in% c("clone", "family","organism", "study","comment", "representative", "short", "type", "filename",
                                                                   "ID", "IC", "length", "consensus","kld_between_revcomp" ))], collapse="_"), "_v2")
      
      
      
    }
    
    PWMs_metadata$filename[m]=filename
    PWMs_metadata$ID[m]=ID
    
    
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE,file=filename, sep="\t")
  
    pcm=as.matrix( read.csv(filename, header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    # Make reverse complement
    # Reverse the matrix = flip the matrix horizontally
    pcm_rev=pcm[,ncol(pcm):1]
    
    #Substitute nucleotides with complements = flip the orders of A&T and C&G
    
    pcm_rev_comp=pcm_rev
    pcm_rev_comp["A", ]=pcm_rev["T",]
    pcm_rev_comp["T", ]=pcm_rev["A",]
    pcm_rev_comp["C", ]=pcm_rev["G",]
    pcm_rev_comp["G", ]=pcm_rev["C",]
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                                profileMatrix=pcm
    )
    
    pfm_rev_comp <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name=PWMs_metadata$symbol[m],
                                         profileMatrix=pcm_rev_comp
    )
    
    pwm_class <- TFBSTools::toPWM(pfm, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    pwm_class_revcomp <- TFBSTools::toPWM(pfm_rev_comp, type="prob", 
                                          pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    # Compute KL divergence, it does not matter in which order you give the matrices to the function
    kl_divergence <- as.numeric(compute_kl_divergence(pwm_class@profileMatrix, pwm_class_revcomp@profileMatrix))
    
    PWMs_metadata[m, "kld_between_revcomp"]=kl_divergence
    
    
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
    PWMs_metadata[m, "IC"]=sum(rowSums(icm)) 
    
    #motif length
    PWMs_metadata[m, "length"]=length(pfm)
  
    motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
    motif@name=ID
    
    PWMs_metadata[m, "consensus"]=motif@consensus
    
    transfac=paste0(pfms_transfac_path,"/", ID,".pfm")
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
  
    write.table(PWMs_list[[m]][,-1],row.names = FALSE, col.names=FALSE, quote=FALSE, 
                file=paste0(pfms_space_path,"/", ID, ".pfm"), sep=" ")
  
    PWM=as.matrix(PWMs_list[[m]][,-1], dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
              append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","/Jolma2015.scpd"))
    append=TRUE
  
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
              file=paste0(pfms_scpd,"/","/Jolma2015.scpd"))
    
    
    #Write also pwms 
    
    write.table(pwm_class@profileMatrix, file=paste0(pwms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pwm_class@profileMatrix, file=paste0(pwms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    #Write reverse complements
    
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
    write.table(pcm_rev_comp, file=paste0(rev_comp_pfms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
    
    
    
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
    
    
    pdf(paste0(logos_pdf_ic,"/",ID,".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    dev.off()
  
  
    #Reverse complements 
    
    png(paste0(revcomp_logos_png_prob,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_prob,"/",ID,".pdf"), width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    png(paste0(revcomp_logos_png_ic,"/",ID,".png"), res=600,  width = 2500/4*ncol(pwm_class_revcomp@profileMatrix), height = 2500)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
    pdf(paste0(revcomp_logos_pdf_ic,"/",ID,".pdf"),  width=ncol(pwm_class_revcomp@profileMatrix), height = 4)
    
    print(ggplot() + ggseqlogo::geom_logo(  pwm_class_revcomp@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
            ggseqlogo::theme_logo()  +
            labs(title=pwm_class_revcomp@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
    
    dev.off()
    
    
  }
}


#Add family mappings
#tmp=PWMs_metadata
TF_family_mappings <- read_delim("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv", delim=",")

#There is no singles in this data
single <- PWMs_metadata$symbol[str_detect(PWMs_metadata$symbol, "_") == FALSE] #31
single_index <- which(PWMs_metadata$symbol %in% single)
Family_names <- rep(NA, nrow(PWMs_metadata))
Family_names[single_index]=as.vector(PWMs_metadata[single_index,"symbol"] %>%
  left_join(select(TF_family_mappings, "symbol", "DBD" ), by = 'symbol') %>% select("DBD"))$DBD


pair <- PWMs_metadata$symbol[str_detect(PWMs_metadata$symbol, "_") == TRUE]
pair_index <- which(PWMs_metadata$symbol %in% pair) #562
#[pair_index]
pairs <- as.data.frame(do.call(rbind, strsplit(PWMs_metadata$symbol[pair_index], split="_")))
colnames(pairs)=c("symbol", "symbol1")
pairs$order=seq(nrow(pairs))

commondf_pairs <- merge(pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,5)]
commondf_pairs <- commondf_pairs[order(commondf_pairs$order), ]
commondf_pairs$order <- NULL

names(commondf_pairs)=c("symbol0", "symbol", "DBD0")
commondf_pairs$order=seq(nrow(commondf_pairs))
commondf_pairs2 = merge(commondf_pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,4,6)]
commondf_pairs2 <- commondf_pairs2[order(commondf_pairs2$order), ]
commondf_pairs2$order <- NULL
commondf_pairs2 <- commondf_pairs2[, c(2,1,3,4)]
names(commondf_pairs2)=c("symbol0", "symbol1", "DBD0", "DBD1")

pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

Family_names[pair_index]=pair_families


PWMs_metadata <- add_column(PWMs_metadata, Lambert2018_families = Family_names, .before = 4)
unique_Lambertfamilies <- unique(PWMs_metadata$Lambert2018_families)

#Is some of the Lamber2018_families missing

which(is.na(PWMs_metadata$Lambert2018_families))

PWMs_metadata[582, 1:5]
PWMs_metadata$Lambert2018_families[582]=PWMs_metadata$family[582]

#Add ensemble gene ID info for monomers







PWMs_metadata <- PWMs_metadata %>%
  left_join(select(TF_family_mappings, "symbol", "Gene Information/ID" ), by = 'symbol')

# PWMs_metadata$`Gene Information/ID`[which(PWMs_metadata$experiment=="HT-SELEX")]

PWMs_metadata <- PWMs_metadata %>% #Save this info only for human monomers
  dplyr::rename("Human_Ensemble_ID"= "Gene Information/ID") %>%
  relocate("Human_Ensemble_ID", .after = "symbol")

#Relocate ID
PWMs_metadata <- PWMs_metadata %>%
  relocate("ID", .before = "symbol")


write.table(PWMs_metadata, file="../../Data/SELEX-motif-collection/Jolma2015_metadata.csv", row.names = FALSE, sep="\t")



