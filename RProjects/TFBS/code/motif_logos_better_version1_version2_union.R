#install.packages("ggseqlogo")
library("ggseqlogo")
library("motifStack")
library("ggplot2")


#make directories

setwd("/projappl/project_2006203/TFBS/RProjects/TFBS")

data_path="/scratch/project_2006203/TFBS"
data_path="/Users/osmalama/projects/TFBS"

dir.create(file.path(data_path, "PWMs_final_union/Logos_final"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final"), "pdf"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final"), "png"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/pdf"), "barcode"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/pdf"), "ic"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/pdf"), "prob"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/png"), "barcode"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/png"), "ic"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final/png"), "prob"), showWarnings = FALSE)

dir.create(file.path(data_path, "PWMs_final_union/Logos_final2"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2"), "pdf"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2"), "png"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/pdf"), "barcode"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/pdf"), "ic"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/pdf"), "prob"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/png"), "barcode"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/png"), "ic"), showWarnings = FALSE)
dir.create(file.path(paste0(data_path, "/PWMs_final_union/Logos_final2/png"), "prob"), showWarnings = FALSE)



pcm_list<-TFBSTools::PFMatrixList()
pwm_list<-TFBSTools::PWMatrixList()


representatives <- read_delim("~/projects/TFBS/PWMs_final_union/metadata_4003_motifs.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)


for(i in 1:nrow(representatives)){
  print(i)
  pcm=as.matrix( read.table(paste0( representatives$filename[i]), header=FALSE) )
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  #A 1294  802 6919 6919    0  156 6919 1801
  #C 6919 3055 1044  249   41 1461  493 2373
  #G 3132  905  938   56    0  536 1133 1625
  #T 2804 3864  657    0 6919 6919  462 1120
  
  
  pcm_class <- TFBSTools::PFMatrix(ID=representatives$ID[i], name=representatives$symbol[i], 
                                   #matrixClass="Zipper-Type", 
                                   strand="+", 
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                   tags=list(family=representatives$Lambert2018_families[i], 
                                             species=representatives$organism[i] 
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
  
  
  
  
  pwm_class <- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  
  
  pcm_list[[representatives$ID[i]]]=pcm_class
  pwm_list[[representatives$ID[i]]]=pwm_class
  
    
  #Probability matrices
  png(paste0(data_path,"/PWMs_final_union/Logos_final2/png/prob/",representatives$ID[i],".png"), res=600,  width = 2500/4*representatives$length[i], height = 2500)
  print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = FALSE, 
                             ylab="probability", font="mono,Courier", fontface = "bold"))
  dev.off()
  
  pdf(paste0(data_path,"/PWMs_final_union/Logos_final2/pdf/prob/",representatives$ID[i],".pdf"), width =representatives$length[i], height = 4, onefile=FALSE)
  print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = FALSE, 
                             ylab="probability", font="mono,Courier", fontface = "bold"))
  dev.off()
  
  png(paste0(data_path,"/PWMs_final_union/Logos_final/png/prob/",representatives$ID[i],".png"), res=600,  width = 2500/4*representatives$length[i], height = 2500)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto"  ) + 
    ggseqlogo::theme_logo()  +
    labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  
  dev.off()
  
  
  pdf(paste0(data_path,"/PWMs_final_union/Logos_final/pdf/prob/",representatives$ID[i],".pdf"), width =representatives$length[i], height = 4, onefile=FALSE)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_medium", col_scheme="auto" ) + 
    ggseqlogo::theme_logo()  +
    labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  dev.off()
  
  #Information content matrices
  
  
  png(paste0(data_path,"/PWMs_final_union/Logos_final2/png/ic/",representatives$ID[i],".png"), res=600,  width = 2500/4*representatives$length[i], height = 2500)
  print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = TRUE, 
                             ylab="bits", font="mono,Courier", fontface = "bold"))
  dev.off()
  
  pdf(paste0(data_path,"/PWMs_final_union/Logos_final2/pdf/ic/",representatives$ID[i],".pdf"), width =representatives$length[i], height = 4, onefile=FALSE)
  print(motifStack::plotMotifLogo( pwm_class@profileMatrix, motifName=pwm_class@name, ic.scale = TRUE, 
                             ylab="bits", font="mono,Courier", fontface = "bold"))
  dev.off()
  
  png(paste0(data_path,"/PWMs_final_union/Logos_final/png/ic/",representatives$ID[i],".png"), res=600,  width = 2500/4*representatives$length[i], height = 2500)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto"  ) + 
    ggseqlogo::theme_logo()  +
    labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  
  dev.off()
  
  
  pdf(paste0(data_path,"/PWMs_final_union/Logos_final/pdf/ic/",representatives$ID[i],".pdf"), width =representatives$length[i], height = 4, onefile = FALSE)
  print(ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits", font="roboto_medium", col_scheme="auto" ) + 
    ggseqlogo::theme_logo()  +
    labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  dev.off()
  
  
  
}

#saveRDS(pcm_list, file=paste0("/scratch/project_2006203/TFBS/RProjects/TFBS/RData/pcm_list_final_union.Rds")
#saveRDS(pwm_list, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/pwm_list_final_union.Rds")


#change the x-axis labels
#plot(motif, xaxis=paste0("pos", seq.int(7)+10))

#plot sequence logo stack
#To show multiple motifs on the same canvas as a sequence logo stack, 
#the distance of motifs need to be calculated first. 
#Previously, MotIV3::motifDistances ( R implementation of STAMP4) is 
#used to calculate the distance. However, 
#The MotIV package were dropped from Bioconductor 3_12. 
#Currently, by default, R implementation of matalign is used. 
#After alignment, users can use plotMotifLogos_finaltack, 
#plotMotifLogos_finaltackWithTree or plotMotifStackWithRadialPhylog 
#to draw sequence Logos_final in different layouts. 
#To make it easy to use, we integrated different functionalities into one workflow function named as motifStack.

#motifStack::motifStack(motifs, layout="stack", ncex=1.0, reorder=FALSE, layout="tree")

#Affinity Logos_final

# plot a sequence logo cloud
# We can also plot a sequence logo cloud for DNA motifs.
# 
# ## assign groups for motifs
# groups <- rep(paste("group",1:5,sep=""), each=10)
# names(groups) <- names(pfms)
# ## assign group colors
# group.col <- brewer.pal(5, "Set3")
# names(group.col)<-paste("group",1:5,sep="")
# ## create a list of pfm objects
# pfms <- mapply(names(pfms), pfms, FUN=function(.ele, .pfm){
#   new("pfm",mat=.pfm, name=.ele)}
#   ,SIMPLIFY = FALSE)
# ## use matalign to calculate the distances of motifs
# hc <- clusterMotifs(pfms)
# ## convert the hclust to phylog object
# library(ade4)
# phylog <- ade4::hclust2phylog(hc)
# ## reorder the pfms by the order of hclust
# leaves <- names(phylog$leaves)
# pfms <- pfms[leaves]
# ## extract the motif signatures
# motifSig <- motifSignature(pfms, phylog, cutoffPval=0.0001, min.freq=1)
# ## draw the motifs with a tag-cloud style.
# motifCloud(motifSig, scale=c(6, .5), 
#            layout="rectangles", 
#            group.col=group.col, 
#            groups=groups, 
#            draw.legend=TRUE)

# representatives$ID[i]


#pdf(paste0("/scratch/project_2006203/TFBS/Logos_final/prob/",name,".pdf"), res=600,  width = 4000, height = 2500)
# png(paste0("../../Logos_final/png/prob/",representatives$ID[i],".png"), res=600,  width = 4000, height = 2500)
# p=ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability" ) + ggseqlogo::theme_logo()  +
#   labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none")
# plot(p)
# dev.off()
# 
# pdf(paste0("../../Logos_final/pdf/prob/",representatives$ID[i],".pdf"), width = 15, height = 5)
# 
# p=ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_slab_bold", col_scheme="auto" ) + ggseqlogo::theme_logo()  +
#   labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none")
# plot(p)
# dev.off
# 
# png(paste0("../../Logos_final/png/prob/",representatives$ID[i],".png"), res=600,  width = 4000, height = 2500)
# p=ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="bits" ) + ggseqlogo::theme_logo()  +
#   labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none")
# plot(p)
# dev.off()
# 
# pdf(paste0("../../Logos_final/pdf/prob/",representatives$ID[i],".pdf"), width = 15, height = 5)
# 
# p=ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix, method="probability", font="roboto_slab_bold", col_scheme="auto" ) + ggseqlogo::theme_logo()  +
#   labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none")
# plot(p)
# dev.off()


# ggseqlogo::list_fonts()
# Available ggseqlogo fonts:
#   helvetica_regular
# helvetica_bold
# helvetica_light
# roboto_medium
# roboto_bold
# roboto_regular
# akrobat_bold
# akrobat_regular
# roboto_slab_bold
# roboto_slab_regular
# roboto_slab_light
# xkcd_regular
# 
# ggseqlogo::list_col_schemes()
# auto
# chemistry
# chemistry2
# hydrophobicity
# nucleotide
# nucleotide2
# base_pairing
# clustalx
# taylor


# font="mono,Courier", fontface="plain"
# 
# png(paste0("../../Logos_final/png/IC/",name,".png"), res=600,  width = 4000, height = 2500)
# p=ggplot() + ggseqlogo::geom_logo(  example_pwms[[name]]@profileMatrix, method="bits" ) + ggseqlogo::theme_logo()  +
#   labs(title=name)+ theme(title =element_text(size=8, face='bold'))
# plot(p)
# #seqLogo(seqLogo_pwms[[name]], fill=c(A="darkgreen", C="blue", G="orange", T="red"))
# dev.off()
# 
# 
# 
# 
# seqLogo::makePWM(pwm_class@profileMatrix)
# 
# install
#ggplot() + ggseqlogo::geom_logo(  pwm_class@profileMatrix ) + ggseqlogo::theme_logo()
