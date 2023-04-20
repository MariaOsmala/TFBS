library(rtracklayer)
#correlation
#install.packages("corrr")
#install.packages("corrplot")
library("corrplot")
library(RColorBrewer)
#install.packages("ClusterR")
library(mclust)
library(tidyr)

representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length.tsv", sep="\t", header=TRUE)
#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
PROXes=representatives[grep("PROX", representatives$symbol),] #9

write.table(PROXes, file="PROXes.csv", sep=",", row.names=FALSE, quote=FALSE)

CAP_selex_PROX_HOX=PROXes[PROXes$experiment=="CAP-SELEX",]
CAP_selex_PROX_HOX=CAP_selex_PROX_HOX[-grep("MAFK", CAP_selex_PROX_HOX$symbol),]

HT_selex_PROX=PROXes[PROXes$experiment!="CAP-SELEX",]


bigbed_filenames=system("ls /scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed*/MOODS_bigbed/*PROX*top.bed", intern=TRUE)

bigbed=gsub("_top.bed", "", do.call(rbind, strsplit(bigbed_filenames, "/"))[,8])
bigbed=gsub("_composite", "", bigbed)
bigbed=gsub("_spacing", "", bigbed)



bigbed_filenames=bigbed_filenames[which( bigbed %in% gsub(".pfm","", do.call(rbind, strsplit( CAP_selex_PROX_HOX$filename, "/"))[,5]) )]


library(R.utils)
sapply(bigbed_filenames,countLines)

#motifs

motifs<-TFBSTools::PFMatrixList()

for(i in 1:nrow(CAP_selex_PROX_HOX)){
  print(i)
  profileMatrix=as.matrix(read.table(file=paste0("/scratch/project_2006203/TFBS/",CAP_selex_PROX_HOX$filename[i])))
  dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
  
  motifs[[CAP_selex_PROX_HOX$ID[i]]] <- TFBSTools::PFMatrix(ID=CAP_selex_PROX_HOX$ID[i], name=CAP_selex_PROX_HOX$symbol[i], 
                                                            #matrixClass="Zipper-Type", 
                                                            strand="+", 
                                                            bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                                            tags=list(family=CAP_selex_PROX_HOX$Lambert2018.families[i], 
                                                                      #species="10090", 
                                                                      #tax_group="vertebrates", 
                                                                      #medline="7592839", 
                                                                      #type="SELEX", 
                                                                      #ACC="P53762", 
                                                                      #pazar_tf_id="TF0000003",
                                                                      T#FBSshape_ID="11", 
                                                                      #TFencyclopedia_ID="580"
                                                            ),
                                                            profileMatrix=profileMatrix
  )
  
  
  
}

PROX_motifs<-TFBSTools::PFMatrixList()

for(i in 1:nrow(HT_selex_PROX)){
  print(i)
  profileMatrix=as.matrix(read.table(file=paste0("/scratch/project_2006203/TFBS/",HT_selex_PROX$filename[i])))
  dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
  
  PROX_motifs[[HT_selex_PROX$ID[i]]] <- TFBSTools::PFMatrix(ID=HT_selex_PROX$ID[i], name=HT_selex_PROX$symbol[i], 
                                                            #matrixClass="Zipper-Type", 
                                                            strand="+", 
                                                            bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                                            tags=list(family=HT_selex_PROX$Lambert2018.families[i], 
                                                                      #species="10090", 
                                                                      #tax_group="vertebrates", 
                                                                      #medline="7592839", 
                                                                      #type="SELEX", 
                                                                      #ACC="P53762", 
                                                                      #pazar_tf_id="TF0000003",
                                                                      T#FBSshape_ID="11", 
                                                                      #TFencyclopedia_ID="580"
                                                            ),
                                                            profileMatrix=profileMatrix
  )
  
  
  
}





example_pwms <- do.call(TFBSTools::PWMatrixList,lapply(PROX_motifs, TFBSTools::toPWM, type="prob",
                                                       pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))

seqLogo_pwms=lapply( lapply(example_pwms, function(x) x@profileMatrix), seqLogo::makePWM)


library("ggplot2")
ggplot() + ggseqlogo::geom_logo(  example_pwms[[1]]@profileMatrix ) + ggseqlogo::theme_logo()


#library(seqLogo)

for(name in names(example_pwms)){
  
  png(paste0("/scratch/project_2006203/TFBS/Logos/IC/",name,".png"), res=600,  width = 4000, height = 2500)
  p=ggplot() + ggseqlogo::geom_logo(  example_pwms[[name]]@profileMatrix, method="bits" ) + ggseqlogo::theme_logo()  +
    labs(title=name)+ theme(title =element_text(size=8, face='bold'))
  plot(p)
  #seqLogo(seqLogo_pwms[[name]], fill=c(A="darkgreen", C="blue", G="orange", T="red"))
  dev.off()
  
  png(paste0("/scratch/project_2006203/TFBS/Logos/prob/",name,".png"), res=600,  width = 4000, height = 2500)
  p=ggplot() + ggseqlogo::geom_logo(  example_pwms[[name]]@profileMatrix, method="probability" ) + ggseqlogo::theme_logo()  +
    labs(title=name)+ theme(title =element_text(size=8, face='bold'))
  plot(p)
  #seqLogo(seqLogo_pwms[[name]], fill=c(A="darkgreen", C="blue", G="orange", T="red"),
  #        ic.scale=FALSE)
  dev.off()
  
}
