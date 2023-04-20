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

MEIS1_DLX3=representatives[grep("MEIS1_DLX3", representatives$filename),]

profileMatrix=as.matrix(read.table(file=paste0("/scratch/project_2006203/TFBS/",MEIS1_DLX3$filename)))
dimnames(profileMatrix)=list(c("A", "C", "G", "T"))


pfm <- TFBSTools::PFMatrix(ID=MEIS1_DLX3$ID, name=MEIS1_DLX3$symbol, 
                                                          #matrixClass="Zipper-Type", 
                                                          strand="+", 
                                                          bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                                          tags=list(family=MEIS1_DLX3$Lambert2018.families, 
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


pwm=TFBSTools::toPWM(pfm, type="prob", pseudocounts=0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
PPMp <- function(C, p){
  (C + p / nrow(C)) / matrix( (colSums(C) + p),nrow=4, ncol=ncol(C),byrow=TRUE)
  
  }

PWM=PPMp(pfm@profileMatrix,0.01)
PWM
pwm@profileMatrix

#Position specific scoring matrix
S <- function(C, B) log2(C / B)

S(PWM, 0.25)

U <- function(C) -colSums(C * log2(C), na.rm = TRUE)
sum(U(PWM)) # 13.01198
sum(2-U(PWM)) # 14.98802 THIS IS THE CORRECT INFORMATION CONTENT

#This is the information content matrix
ICM=(2-matrix( U(PWM),nrow=4, ncol=ncol(C),byrow=TRUE))*PWM


icm=TFBSTools::toICM(pfm, pseudocounts=0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))


TFBSTools::seqLogo(icm)

MEIS1_DLX3

setwd("/home/local/osmalama/projects/TFBS/RProjects/TFBS/")
PCM <- importMatrix("../../PWMs/Jolma2015/pwms/Homo_sapiens/MEIS1_DLX3_CAP-SELEX_TAAAGC40NGAA_AT_TGACANSNTAATTG_1_3_NO.pfm")




#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
PROXes=representatives[grep("PROX", representatives$symbol),] #9

CAP_selex_PROX_HOX=PROXes[PROXes$experiment=="CAP-SELEX",]
CAP_selex_PROX_HOX=CAP_selex_PROX_HOX[-grep("MAFK", CAP_selex_PROX_HOX$symbol),]

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

example_pwms <- do.call(TFBSTools::PWMatrixList,lapply(motifs, TFBSTools::toPWM, type="prob",
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



#The first option controls colors used for plotting BIC, ICL, etc. curves, 
#whereas the second option is used to assign colors for indicating clusters or classes when plotting data.
library("mclust")
mclust.options("bicPlotColors")
mclust.options("classPlotColors")

#Starting with R version 4.0, the function can be used for retrieving colors from some pre-defined palettes. For instance
cbPalette=palette.colors(palette = "Okabe-Ito")
bicPlotColors <- mclust.options("bicPlotColors")
bicPlotColors[1:14] <- c(cbPalette, cbPalette[1:5])
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = cbPalette[-1])


#TF bedfile
#TF=/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/CTCF_HT-SELEX_TAGCGA20NGCT_AJ_NGCGCCMYCTAGYGGTN_2_4_YES_top.bed
#wc -l $TF #51136
#TF=import("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_3_processed/MOODS_bigbed/CTCF_HT-SELEX_TAGCGA20NGCT_AJ_NGCGCCMYCTAGYGGTN_2_4_YES_top.bed")
#77283
#117589
#phylop=/scratch/project_2007567/phyloP/241-mammalian-2020v2.bigWig

#PROX1_HOXA2



TF=import(bigbed_filenames[grep("PROX1_HOXA2",bigbed_filenames)])

#One of the most common ways of storing score data is, as mentioned, the wig or bigWig format. 
#In addition, conservation scores can also be downloaded in the wig/bigWig format. 
#You can import bigWig files into R using the import() function from the rtracklayer package. 
#However, it is generally not advisable to read the whole bigWig file in memory 
#as was the case with BAM files. Usually, you will be interested in only a 
#fraction of the genome, such as promoters, exons etc. S
#o it is best that you extract the data for those regions and 
#read those into memory rather than the whole file. 
#Below we read a bigWig file only for the bases on promoters. 
#The operation returns a GRanges object with the score column which indicates the scores in the bigWig file per genomic region.


phylop.gr=import("/scratch/project_2007567/phyloP/241-mammalian-2020v2.bigWig", which=TF,as = "RleList") # get coverage vectors

#6.3.1 Extracting subsections of Rle and RleList objects
#Frequently, we will need to extract subsections of the Rle vectors or RleList objects. 
#We will need to do this to visualize that subsection or get some statistics out of those sections. 
#For example, we could be interested in average coverage per base 
#for the regions we are interested in. 
#We have to extract those regions from the RleList object and apply summary statistics. 
#Below, we show how to extract subsections of the RleList object.
#We are extracting promoter regions from the ChIP-seq read coverage RleList. 
#Following that, we will plot one of the promoter’s coverage values.


myViews=Views(phylop.gr,as(TF,"IRangesList")) # get subsets of coverage
# there is a views object for each chromosome
myViews
myViews[[1]]

plot(myViews[[1]][[5]],type="l") #This is just one row in a single chromosome

head(
  viewMeans(myViews[[1]])
)


phyloP=GRanges(myViews)
phyloP$view=NULL

sm=ScoreMatrix(target=phyloP, windows=TF, strand.aware=TRUE, weight.col="score")

#could be computed also as a rowMeans of sm
rwmns=rowMeans(sm@.Data)


#Figure out the consensus sequence for each TF, the most probable sequence
df=data.frame(cm=colMeans(sm@.Data))
ggplot(df, aes(y=cm, x=1:17))+geom_line()+xlab("")+ylab("Average PhyloP score")+ theme_minimal()+labs(title="CTCF")
  scale_x_continuous(name ="", breaks=1:17, labels=c("A","G", "C", "G", "C", "C","A", "C", "C", "T", "A", "G", "T", "G", "G", "T", "N"))+
  theme(panel.grid.minor = element_blank())
                   
#heatMatrix(sm, xcoords = seq(1:ncol(sm) ))
#plotMeta(sm, xcoords  = seq(1:ncol(sm) ), ylab = "average PhyloP score")
#plotMeta(sm, xcoords  = c("A","G", "C", "G", "C", "C","A", "C", "C", "T", "A", "G", "T", "G", "G", "T", "N"), ylab = "average PhyloP score")


# Using different color spectrum
M=cor(sm@.Data)
col<- colorRampPalette(c("red", "white", "blue"))(20)
corrplot(M, method="shade", type="lower", diag=FALSE, col=col,tl.col="black", tl.srt=45) #"circle" "square" "ellipse" "number" "shade" "color" "pie"


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="shade", type="lower", diag=FALSE, col=brewer.pal(n=10, name="RdBu"),addCoef.col = "black") # Add coefficient of correlation) 

corrplot(M, method="shade", type="lower", diag=FALSE, col=col(200),addCoef.col = "black")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat<-cor.mtest(sm@.Data)

corrplot(M, type="upper", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

#bg="lightblue" #background color if you use circles or squares

#If x and y are matrices then the covariances (or correlations) between the columns of x and the columns of y are computed.


  colpair_map(.data, .f, ..., .diagonal = NA)
  
  
conservation_means=unlist(viewMeans(myViews))

#correlations between conservations

#all pair-wise correlations between motif positions
str(sm@.Data)


#nstart	if centers is a number, how many random sets should be chosen?

#km_result <- kmeans(conservation_means, centers = 2, iter.max = 10, nstart = 20)
#The following models are available in package mclust: univariate mixture 
# "E" equal variance (one-dimensional) 
# "V" variable/unqual variance (one-dimensional)

BIC <- mclustBIC(conservation_means)
plot(BIC)

gm_result_2 <- Mclust(conservation_means, G=2, modelNames="V")
gm_result_3 <- Mclust(conservation_means, G=3, modelNames="V")
gm_result_4 <- Mclust(conservation_means, G=4, modelNames="V")
#gmda_result <- MclustDA(conservation_means, gm_result$classification, modelType="EDDA")

fit_2 <- densityMclust(conservation_means, G=2)
fit_3 <- densityMclust(conservation_means, G=3)

summary(gm_result, parameters = TRUE)
table( gm_result$classification)

summary(fit)
plot(fit, what = "density", data = conservation_means, breaks = 100)
plot(fit, what = "diagnostic", type = "cdf")
plot(fit, what = "diagnostic", type = "qq")
plot(gm_result, what = "classification", main = "Mclust Classification")
plot(gm_result, what = "uncertainty")
plot(gm_result, what = "density")
#plot(gmda_result, what = "scatterplot")
#summary(gmda_result)

df.pp<-data.frame(PhyloP = conservation_means)

ggplot(df.pp, aes(x = PhyloP)) +
  geom_histogram(binwidth = 0.05, alpha=0.5) +
  mapply(
    function(mean, sd, lambda, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
        }
      )
    },
    mean = gm_result_2$parameters$mean, #mean
    sd = sqrt(gm_result_2$parameters$variance$sigmasq), #standard deviation
    n = length(df.pp$PhyloP), #sample size
    lambda=gm_result_2$parameters$pro,
    binwidth = 0.05 #binwidth used for histogram
  )

conserved=which(gm_result_2$classification==which.max(gm_result$parameters$mean)) #33306
nonconserved=which(gm_result_2$classification== which.min(gm_result$parameters$mean)) #84283

#PhyloP scores for conserved


