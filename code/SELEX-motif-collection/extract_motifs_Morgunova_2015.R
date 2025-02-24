

#############Read excel tables################

library(tidyverse)
library(readxl)
library(universalmotif)
library(ggplot2)
rm(list=ls())
source("read_excel_tables.R")
source("paths_to_folders.R", echo=FALSE)
source("compute_kl_divergence.R", echo=FALSE)
#Adjacent dinucleotide model
#17. Jolma, A. et al. DNA-binding specificities of human transcription factors. Cell
#152, 327–339 (2013).

dinucleotide=read.table("../../Data/SELEX-motif-collection/E2F8_cycle4.adm", header=FALSE,nrows = 16)

families=read_csv("../../Data/SELEX-motif-collection/HumanTFs-ccbr-TableS1.csv")
#tmp=families[grep("E2F8", families$symbol),"DBD"]$DBD


pwm=read.table("../../Data/SELEX-motif-collection/E2F8_cycle4.adm", header=FALSE, skip = 16)
pwm=pwm[,-c(15,16)]

pcm=as.matrix(pwm)*1000
dimnames(pcm)=list(c("A", "C", "G", "T"))

# Make reverse complement
# Reverse the matrix = flip the matrix horizontally
pcm_rev=pcm[,ncol(pcm):1]

#Substitute nucleotides with complements = flip the orders of A&T and C&G

pcm_rev_comp=pcm_rev
pcm_rev_comp["A", ]=pcm_rev["T",]
pcm_rev_comp["T", ]=pcm_rev["A",]
pcm_rev_comp["C", ]=pcm_rev["G",]
pcm_rev_comp["G", ]=pcm_rev["C",]


pcm_class <- TFBSTools::PFMatrix(ID="E2F8", name="E2F8", 
                                strand="+", 
                                 bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                 tags=list(family=families[grep("E2F8", families$symbol),"DBD"]$DBD, 
                                           species="Homo_sapiens" 
                                 ),
                                 profileMatrix=pcm
)

pcm_rev_comp <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),name="E2F8",
                                     profileMatrix=pcm_rev_comp
)


pwm_class <- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

pwm_class_revcomp <- TFBSTools::toPWM(pcm_rev_comp, type="prob", 
                                      pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

#Compute the KL-divergence between the original and the reverse complement

# Compute KL divergence, it does not matter in which order you give the matrices to the function
kl_divergence <- as.numeric(compute_kl_divergence(pwm_class@profileMatrix, pwm_class_revcomp@profileMatrix))




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







# Draw logos --------------------------------------------------------------



png(paste0(logos_png_prob,"/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))

dev.off()


pdf(paste0(logos_pdf_prob,"/","E2F8_Morgunova2015",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
dev.off()

#Information content matrices


png(paste0(logos_png_ic,"/","E2F8_Morgunova2015",".png"), res=600,  width = 2500/4*ncol(pwm_class@profileMatrix), height = 2500)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))

dev.off()


pdf(paste0(logos_pdf_ic,"/","E2F8_Morgunova2015",".pdf"), width=ncol(pwm_class@profileMatrix), height = 4)
print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
        ggseqlogo::theme_logo()  +
        labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
dev.off()

ID="E2F8_Morgunova2015"

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



#pfms_tab_path
#pfms_space_path
#pfms_transfac_path
#pfms_scpd
#pwms_tab_path
#pwms_space_path

write.table(pwm_class@profileMatrix, file=paste0(pwms_space_path,"/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
write.table(pwm_class@profileMatrix, file=paste0(pwms_tab_path,"/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 

write.table(pcm_class@profileMatrix, file=paste0(pfms_tab_path,"/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
write.table(pcm_class@profileMatrix, file=paste0(pfms_space_path,"/","E2F8_Morgunova2015", ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 

#Write reverse complements

write.table(pcm_rev_comp@profileMatrix, file=paste0(rev_comp_pfms_space_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
write.table(pcm_rev_comp@profileMatrix, file=paste0(rev_comp_pfms_tab_path,"/",ID, ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 


icm <- TFBSTools::toICM(pcm_class, pseudocounts=0.01, schneider=FALSE)
sum(TFBSTools::rowSums(icm))

#motif length
length(pcm_class)

filename=paste0(pfms_tab_path,"/","E2F8_Morgunova2015", ".pfm")
motif=universalmotif::read_matrix(file=filename, sep="\t", header=FALSE)
motif@name="E2F8_Morgunova2015"

tmp=data.frame(ID="E2F8_Morgunova2015",
               symbol="E2F8",
               Human_Ensemble_ID=families[grep("E2F8", families$symbol),]$`Gene Information/ID`, 
               clone="", 
               family="E2F",
               Lambert2018_families="E2F",
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
               comment=families[grep("E2F8", families$symbol),]$`Binding mode`,
               filename=filename,
               IC=sum(TFBSTools::rowSums(icm)),
               length=length(pcm_class),
               consensus=motif@consensus, 
               kld_between_revcomp=kl_divergence

)

write.table(tmp, file="../../Data/SELEX-motif-collection/Morgunova2015_metadata.csv", row.names = FALSE, sep="\t")


transfac=paste0(pfms_transfac_path, "/", tmp$ID,".pfm")
write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)


PWM=as.matrix(pcm, dimnames=NULL)
rownames(PWM)=c("A", "C", "G", "T")
write.table(paste0(">",  tmp$ID),   
            append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE,
            file=paste0(pfms_scpd,"/","/Morgunova2015", ".scpd"))
append=TRUE

write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
            file=paste0(pfms_scpd,"/","/Morgunova2015", ".scpd"))






