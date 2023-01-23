library("rtracklayer")
library("dbplyr")
library("dplyr")

#short names for TFs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))

which(representatives$ID %in% rownames(p_matrix))

shorter_names=representatives %>% filter(test %in% rownames(p_matrix)) %>% select(symbol)
table(shorter_names)

#These are log_10(p-values)
p_matrix=readRDS(  file = paste0( "../RData/Vierstra_p_matrix_human_cell_type_restricted.Rds")) #

#FDR correction for the p-values?

#These are log2( (k / n) / (m / N) )
e_matrix=readRDS(  file = paste0( "../RData/Vierstra_e_matrix_human_cell_type_restricted.Rds")) #

#Heatmaps showing Gene Ontology TF motifs with maximal enrichment in cell-type-restricted cCREs of selected cell types. 
#Only the most enriched TF motif in each of the previously identified motif archetypes (Vierstra et al., 2020) was selected as the representative and 
#the top 10 motifs were selected for each cell type. Color represents −log10P.

#which are enriched 

infinites=which(is.infinite(e_matrix))

e_matrix[infinites]=NA

enriched=which(e_matrix>0) #1062*111=117882
#50228

#why this does not work
# r = ((enriched-1) %% ncol(e_matrix)) + 1
# c = floor((enriched-1) / ncol(e_matrix)) + 1
# e_matrix[r,c]

enrichments=e_matrix[enriched] #-Inf to 2.7
min(enrichments, na.rm=TRUE) #1.417561e-05
max(enrichments, na.rm=TRUE) #2.711503

#which are depleted 
depleted=which(e_matrix<0)
depletions=e_matrix[depleted]
min(depletions, na.rm=TRUE) #-2.734187
max(depletions, na.rm=TRUE) #-8.023386e-06

#exactly zero?
#equal=which(e_matrix==0) empty


#Convert the depleted -log_10 p-values to negative
str(depleted)
p_matrix_depleted=-p_matrix
p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]

#plot heatmap
library("pheatmap")
library("cluster")

#cluster the heatmap

#dim(p_matrix_depleted)
ed <- dist(t(p_matrix_depleted)) #Euclidean distance is the default, N times D matrix as input

# do the hierarchical clusterings
#method: the agglomeration method to be used. 
#This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#"single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

#This clusters now the cell lines based on TF binding
ehc <- hclust(ed,method="complete")
ehc2 <- hclust(ed,method="single")
ehc3 <- hclust(ed,method="ward.D")
ehc4 <- hclust(ed,method="ward.D2")
ehc5 <- hclust(ed,method="average")
ehc6 <- hclust(ed,method="mcquitty")
ehc7 <- hclust(ed,method="median")
ehc8 <- hclust(ed,method="centroid")
# plot the resulting dendrogram: complete linkage method produces clear groups (2 or 3 clear clusters), 
#while single linkage starts from a single patient and goes to the most dissimilar which is not very informative here.
plot(ehc,labels=FALSE)
plot(ehc2,labels=FALSE)
plot(ehc3,labels=FALSE)#this looks nice, ward.D
plot(ehc4,labels=FALSE)
plot(ehc5,labels=FALSE)
plot(ehc6,labels=FALSE)
plot(ehc7,labels=FALSE)
plot(ehc8,labels=FALSE)

## Silhouette measure plots
#Display the silhouette measure for the Clustering result based on Euclidean distance 
#for $K=2,\cdots,5$ clusters on the same plot. What is the optimal number of cluster based on the silhouette plots? 
#Hint: Search online on the meaning of silhouette plots.

#Silhouette refers to a method of interpretation and validation of consistency within clusters of data. 
#The technique provides a succinct graphical representation of how well each object has been classified.[1] 

#The silhouette value is a measure of how similar an object is to its own cluster (cohesion) 
#compared to other clusters (separation). The silhouette ranges from −1 to +1, 
#where a high value indicates that the object is well matched to its own cluster 
#and poorly matched to neighboring clusters. If most objects have a high value, 
#then the clustering configuration is appropriate. If many points have a 
#low or negative value, then the clustering configuration may have too many or too few clusters.



#We will do three things in one nested command:
#First (inner-most brackets), find the clusters using cutree(). 
#Second, calculate the silhouette measure with silhouette(). Third (outer-most brackets), make the plot.


# to demonstrate the principle, we can run the last command (for 5 clusters) in three separate commands (unhashtag the lines to see how it works)
ct <- cutree(ehc3,5)
sh <- silhouette(ct,ed)
plot(sh)

plot(silhouette(cutree(ehc3,2),ed))
plot(silhouette(cutree(ehc3,3),ed))
plot(silhouette(cutree(ehc3,4),ed))
plot(silhouette(cutree(ehc3,5),ed))
plot(silhouette(cutree(ehc3,6),ed))
plot(silhouette(cutree(ehc3,7),ed))
plot(silhouette(cutree(ehc3,8),ed))
plot(silhouette(cutree(ehc3,9),ed))
plot(silhouette(cutree(ehc3,10),ed))


#Display the results by a heatmap. Use the R function pheatmap for your subset of 1,500 genes with the most variance from Task 2.1. 
#Add annotations of the first four mutations found in the clinical data. 
#Use "complete" as the clustering method. Include the optimal number of clusters found in Task 2.3. 
#What is visualized in addition to the dendrogram on your heatmap? 
#  ```{r heatmap, fig.cap = 'Pretty Heatmap with Sample Annotation', warning=FALSE}
# generate an annotation data.frame to visualize the first four mutations
#anno <- data.frame(clinic[,27:30]) #176 x 90

# plot heatmap with sample annotation. Heat map produces both clustering info and gene expression info. The colors correspond to the gene expression levels. Additional annotation allows to see certain genes' mutation status in patients.

min(p_matrix_depleted)
max(p_matrix_depleted)

library("RColorBrewer")
paletteLength=100
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdBu")))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))

pdf(file = "../Figures/adult-cCREs-heatmap.pdf", width=7, height=15 )
pheatmap(t(p_matrix_depleted), color=color,clustering_method="ward.D", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE,
         show_colnames=FALSE,legend=TRUE, breaks= myBreaks, cutree_rows=6, 
          )
dev.off()

#cutree_cols=2)

#Epithelial
#Visualise the Hepatocyte, Acinar, Parietal, Chief, Foveolar, G. Neuroendo, Ductal, 

rownames(p_matrix_depleted)=shorter_names$symbol

epithelial=p_matrix_depleted[, which( colnames(p_matrix_depleted) %in% 
                                        c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal") )]
       
#what are the 10 most enriched motifs for these cell lines       

most_enriched<-list()
for(c in colnames(epithelial)){
  
  most_enriched[[c]]=names(sort(epithelial[,c], decreasing = TRUE)[1:10])
  
}
 
#41  
most_enriched=unique(unlist(most_enriched))
   
epithelial=epithelial[most_enriched,] 


epithelial[
  order( epithelial[,1], epithelial[,2],epithelial[,3],epithelial[,4],epithelial[,5],epithelial[,6] , decreasing=FALSE),
]

epithelial[
  with(epithelial, order("Foveolar", decreasing=FALSE)),
]

epithelial=epithelial[
  order( epithelial[,c( "Foveolar")] ),
]

epithelial=epithelial[,c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal")]

epithelial=epithelial[order(rownames(epithelial)),]

epithelial=epithelial[c("HNF1A","HNF4A","ESRRA","GATA2","GATA3","GATA5", "FOXA1", "FOXA2" , "FOXC2",  "FOXD2",  "ONECUT1",  "ONECUT2",  "NR2F1", "BARHL1_TEAD4", "CEBPB_FOXD2",
  "CEBPG",  "FIGLA",  "FOSL1", "GSC2_TEAD4", "HOXA10_TEAD4", "HOXB2_TCF3", "HOXD12_TEAD4", "ID4", "JUN", "LMX1A_BACH2",  "MYBL1_FIGLA",  "NEUROG1_TCF3",
"NR2C1", "OTX1_TEAD4", "POU4F1", "RARA","RFX1", "RFX3",  "RFX3_FIGLA", "TCF12", "TEAD2", "TFAP4_FLI1" ),]

library(grid)




source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")


pdf(file = "../Figures/adult-epithelial-heatmap.pdf", width=7, height=4 )

pheatmap(t(epithelial), color=color, clustering_method="ward.D",legend=TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks= myBreaks, row_names_side = "left", row_dend_side = "right")
dev.off()
       
# set up panel of 2x2 plots:
par(mfrow=c(2,2))

# plots for 2, 3, and 4 clusters:



#Which are representative motifs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]
rep_motifs2=rep_motifs

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))
#1062



cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs)


N=length(all_cCREs) #all possible cCREs #435142



library(ggplot2)
# Basic histogram

df=data.frame(cCREs_nro= unlist(lapply(cCREs_list, length)))

max(df$cCREs_nro) #48391
min(df$cCREs_nro) #846

#ggplot(df, aes(x=_nro)) + geom_histogram()
# Change the width of bins
pdf(file = "../Figures/cCREs-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=cCREs_nro)) + geom_histogram(binwidth=1000)
dev.off()

#Width of cCREs, these are all the same 400 bp?

widths=lapply(cCREs_list, width)
lapply(widths, table)

match_numbers=lapply(representative_motif_matches_Vierstra, length)

df2=data.frame(match_nro=unlist(match_numbers))

max(df2$match_nro) #4063189
min(df2$match_nro) #21108


# Change the width of bins
pdf(file = "../Figures/Vierstra-representative-motif-hit-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df2, aes(x=match_nro)) + geom_histogram(binwidth=10000)
dev.off()

mean(df2$match_nro) #307578
median(df2$match_nro) #247368

#Are the low hitters longer motifs/dimers and super hitters short

#select only 1062 representatives
representatives_match_numbers=representatives[which(representatives$new_representative=="YES"),]
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives_match_numbers$filename,"/"))[,5])

table(test==rownames(df2) ) #the motif names are in the same order

representatives_match_numbers$match_numbers=df2$match_nro

monomers=c("monomeric", "monomer or dimer", "monomer", "monomeric or dimeric")

multimers=c("dimeric", "putative multimer", "putatively multimeric", "dimer of dimers", "trimeric")


#both monomeric

representatives_match_numbers$motif_type=rep("monomer", nrow(representatives_match_numbers))

cap_selex=representatives_match_numbers$experiment=="CAP-SELEX"

mm=representatives_match_numbers$type %in% multimers

representatives_match_numbers[mm | cap_selex, "motif_type"]="multimer"

# Overlaid histograms
ggplot(representatives_match_numbers, aes(x=match_numbers, color=motif_type)) + geom_histogram(fill="white",binwidth=10000)+ xlim(0, 1e6)

#There seems to be not much difference between the distributions between monomers and multimers

#Means and medians of monomers and multimers

representatives_match_numbers %>%
  group_by(motif_type) %>%
  dplyr::summarize(Mean = mean(match_numbers, na.rm=TRUE), Median = median(match_numbers, na.rm=TRUE))

#lengths of the motifs

representatives_match_numbers$motif_lengths=lapply(representative_motif_matches_Vierstra, function(x) unique(width(x)) )

#count as the function of motif length

ggplot(representatives_match_numbers, aes(x=as.numeric(motif_lengths), y=as.numeric(match_numbers))) + geom_point()+ geom_smooth()

+scale_x_continuous()+scale_y_continuous()


#how many have more than half million

length(which(representatives_match_numbers$match_numbers>500000)) #128


#what are these motifs
super_hitters=representatives_match_numbers[representatives_match_numbers$match_numbers> 500000,]
  
  
#What are those with low number of hits
low_hitters=representatives_match_numbers[representatives_match_numbers$match_numbers< 50000,]


#distribution of scores for high hitters and low hitters

#the highest hitter
arrange(super_hitters, desc(match_numbers))$ID[1]
#why there is motif ID, remove
hist(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score)

#How many if using score threshold 5
table(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score>=5)

hist(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score[representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score>=5])

#the lowest hitter
arrange(low_hitters, match_numbers)$ID[1]
hist(representative_motif_matches_Vierstra[[ arrange(low_hitters, match_numbers)$ID[1] ]]$score)


for(hits in motifs){
  print(hits)
  #hits=motifs[1]

 
  #m: the size of all cCREs that have a motif match#
  fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  
  m=length(unique(fol@to)) #This is now number of all cCREs that have at least one motif hit #9857
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #findOverlaps(query, subject)
    fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
    #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
    # number of cCREs that have motif match
    
    k=length(unique(fol@to)) #This is now number of cCREs that have at least one motif hit #178
    #k=length(unique(fol@from)) This is now number of hits that overlap with cCREs (there can be multiple motifs hits in a single CRE)
    
    #m	the number of white balls in the urn.
    # number of motif matches, 300 000
    
    # N-m
    # the number of black balls in the urn.
    # the number of all possible cCREs - the number of motif matches
    
    #Nm=N-m #list of all CREs that do not have a motif This is computed in FishersExact
    
    # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
    # The number of cell type specific cCREs
    
    n=length(cCREs_list[[cell_type]]) #This is correct
    
    #Hypergeometric test/Fisher's exact test
    
    #m: a list of all cCREs that have a motif match
    
    #n: A set of cCREs restrictive to a given cell type. G_0, n=|G_0|
    
    #N: G: all cCREs across all cell types G, N=|G|
    
    #k: number of cell-line-specific cCREs that have motif match
    
    
    fe=FishersExact(k,n,m,N)
    p_matrix[hits, cell_type]=fe$p
    e_matrix[hits, cell_type]=fe$e
    
    
    
  }
  
}


saveRDS(p_matrix,  file = paste0( "../RData/Vierstra_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/Vierstra_e_matrix_human_cell_type_restricted.Rds")) #





