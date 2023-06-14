

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
rm(list=ls())

#Adjacent dinucleotide model
#17. Jolma, A. et al. DNA-binding specificities of human transcription factors. Cell
#152, 327–339 (2013).

dinucleotide=read.table("../../PWMs/Morgunova2015/E2F8_cycle4.adm", header=FALSE,nrows = 16)

families=read_csv("../../../motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv")
tmp=families[grep("E2F8", families$symbol),"DBD"]$DBD

unique(families$DBD)
pwm=read.table("../../PWMs/Morgunova2015/E2F8_cycle4.adm", header=FALSE, skip = 16)
pwm=pwm[,-c(15,16)]

pcm=as.matrix(pwm)*1000
dimnames(pcm)=list(c("A", "C", "G", "T"))

pcm_class <- TFBSTools::PFMatrix(ID="E2F8", name="E2F8", 
                                 #matrixClass="Zipper-Type", 
                                 strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                 tags=list(family=families[grep("E2F8", families$symbol),"DBD"]$DBD, 
                                           species="Homo_sapiens" 
                                           #tax_group="vertebrates", 
                                           #medline="7592839", 
                                           #type="SELEX", 
                                           #ACC="P53762", 
                                           #pazar_tf_id="TF0000003",
                                           #TFBSshape_ID="11", 
                                           #TFencyclopedia_ID="580"
                                 ),
                                 profileMatrix=pcm
)




pwm_class <- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


#The two matrices below are the same
pwm_class@profileMatrix
#A 0.314 0.249 0.305 0.160 0.028 0.002 0.004 0.001 0.001 0.003 0.774 0.676 0.446 0.244
#C 0.121 0.127 0.094 0.052 0.110 0.007 0.989 0.001 0.219 0.033 0.066 0.048 0.157 0.185
#G 0.231 0.213 0.176 0.231 0.850 0.989 0.003 0.996 0.776 0.938 0.060 0.170 0.159 0.258
#T 0.334 0.411 0.425 0.557 0.012 0.002 0.004 0.002 0.004 0.026 0.100 0.106 0.238 0.313
#pwm
#1 0.314 0.249 0.305 0.160 0.028 0.002 0.004 0.001 0.001 0.003 0.775 0.677 0.446 0.244
#2 0.121 0.127 0.094 0.052 0.110 0.007 0.988 0.001 0.219 0.033 0.066 0.048 0.157 0.185
#3 0.231 0.213 0.176 0.231 0.850 0.990 0.003 0.996 0.776 0.938 0.060 0.170 0.159 0.258
#4 0.334 0.411 0.425 0.557 0.012 0.002 0.004 0.002 0.004 0.026 0.100 0.106 0.238 0.313



#Probability matrices
png(paste0("../../Logos2/png/prob/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = FALSE, 
                                 ylab="probability", font="mono,Courier", fontface = "bold"))
dev.off()

pdf(paste0("../../Logos2/pdf/prob/","E2F8_Morgunova2015",".pdf"), width =ncol(pwm_class@profileMatrix), height = 4)
print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = FALSE, 
                                 ylab="probability", font="mono,Courier", fontface = "bold"))
dev.off()

library("ggplot2")
png(paste0("../../Logos/png/prob/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))

dev.off()


pdf(paste0("../../Logos/pdf/prob/","E2F8_Morgunova2015",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
dev.off()

#Information content matrices


png(paste0("../../Logos2/png/ic/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = TRUE, 
                                 ylab="bits", font="mono,Courier", fontface = "bold"))
dev.off()

pdf(paste0("../../Logos2/pdf/ic/","E2F8_Morgunova2015",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = TRUE, 
                                 ylab="bits", font="mono,Courier", fontface = "bold"))
dev.off()

png(paste0("../../Logos/png/ic/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))

dev.off()


pdf(paste0("../../Logos/pdf/ic/","E2F8_Morgunova2015",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
dev.off()

#Write pcm and pfm space and tab separated
pwms="../../PWMs_final/Morgunova2015/pwms/Homo_sapiens/"
pwms_space="../../PWMs_final/Morgunova2015/pwms_space/Homo_sapiens/"
pwms_transfac="../../PWMs_final/Morgunova2015/pwms_transfac/Homo_sapiens/"
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(pwms_transfac, recursive=TRUE)

write.table(pwm_class@profileMatrix, file=paste0("../../PFMs_space/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
write.table(pwm_class@profileMatrix, file=paste0("../../PFMs_tab/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 

write.table(pcm_class@profileMatrix, file=paste0("../../PWMs_final/Morgunova2015/pwms/Homo_sapiens/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
write.table(pcm_class@profileMatrix, file=paste0("../../PWMs_final/Morgunova2015/pwms_space/Homo_sapiens/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 





icm <- TFBSTools::toICM(pcm_class, pseudocounts=0.01, schneider=FALSE)
sum(TFBSTools::rowSums(icm)) #6.77892

#motif length
length(pcm_class)

filename=paste0("../../PWMs_final/Morgunova2015/pwms/Homo_sapiens/","E2F8_Morgunova2015", ".pfm")
motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
motif@name="E2F8_Morgunova2015"

PWMs_metadata[m, "IC_universal"]=motif@icscore
PWMs_metadata[m, "consensus"]=motif@consensus


representatives <- read_tsv("../../../motif-clustering-Viestra-private/metadata/new_representatives_IC_length.tsv", col_names=TRUE)

tmp=data.frame(ID="E2F8_Morgunova2015",
               symbol="E2F8",
               clone="", 
               family="E2F",
               Lambert2018.families="E2F",
               organism="Homo_sapiens",            
               study="Morgunova2015", 
               experiment="HT-SELEX",
               ligand="",
               batch="", 
               seed="",  
               multinomial="",         
               cycle="", 
               representative=NA,
               short="", 
               type="",  
               comment="",
               filename="../../PWMs_final/Morgunova2015/pwms/Homo_sapiens/E2F8_Morgunova2015.pfm",
               IC=sum(TFBSTools::rowSums(icm)),
               length=length(pcm_class),
               IC_universal=motif@icscore,
               consensus=motif@consensus

)

write.table(tmp, file="../../PWMs_final/Morgunova2015/metadata.csv", row.names = FALSE, sep="\t")
saveRDS(tmp, file="Rdata/Morgunova2015.Rds")

transfac=paste0(pwms_transfac, "/", tmp$ID,".pfm")
write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)


PWM=as.matrix(pcm, dimnames=NULL)
rownames(PWM)=c("A", "C", "G", "T")
write.table(paste0(">",  tmp$ID),   
            append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE,
            file=paste0("../../PWMs_final/Morgunova2015/all", ".scpd"))
append=TRUE

write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
            file=paste0("../../PWMs_final/Morgunova2015/all", ".scpd"))





#representatives
  
representatives=bind_rows(representatives, as_tibble(tmp))

write_tsv(representatives, file="../../../motif-clustering-Viestra-private/metadata/new_representatives_IC_length_Morgunova.tsv")


