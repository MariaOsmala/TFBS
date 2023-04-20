
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("motifStack")
BiocManager::install("TFBSTools")

library("motifStack")


MEIS1_DLX3

RUNX1 <- importMatrix(system.file("extdata", "MA0002.1.jaspar",
                                  package = "motifStack",
                                  mustWork = TRUE))[[1]]
setwd("/home/local/osmalama/projects/TFBS/RProjects/TFBS/")
PCM <- importMatrix("../../PWMs/Jolma2015/pwms/Homo_sapiens/MEIS1_DLX3_CAP-SELEX_TAAAGC40NGAA_AT_TGACANSNTAATTG_1_3_NO.pfm")

PFM=(pcm2pfm(PCM)[[1]])

sum(motifStack::getIC(PCM))

#score for sequence

sequence=c("T","G","A","C","A","G","C","A","T","A","A","T","T","G")
sum(diag(log(PFM@mat/0.25)[sequence,]))

134+709+44+618+10+865+929+705

sum(log(PFM@mat[sequence,1:14]/0.25))

plot(PWM[[1]])


library(seqLogo)

p <- makePWM(PFM@mat)

png("../../Figures/ic-logo.png", res=600,  width = 4000, height = 2500)
seqLogo(p, fill=c(A="darkgreen", C="blue", G="orange", T="red"))
dev.off()

png("../../Figures/prob-logo.png", res=600,  width = 4000, height = 2500)
seqLogo(p, fill=c(A="darkgreen", C="blue", G="orange", T="red"),
        ic.scale=FALSE)
dev.off()

#Information content
pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
sum(getIC(motif))

