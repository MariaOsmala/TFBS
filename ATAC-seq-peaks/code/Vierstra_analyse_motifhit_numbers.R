library("rtracklayer")
library("dbplyr")
library("dplyr")



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

representative_motif_matches_Vierstra <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_Vierstra.Rds")

motifs=names(representative_motif_matches_Vierstra)


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
match_numbers=readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Vierstra.Rds")


df2=data.frame(match_nro=unlist(match_numbers))

max(df2$match_nro) #4063189
min(df2$match_nro) #21108


# Change the width of bins
pdf(file = "/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/Vierstra-representative-motif-hit-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df2, aes(x=match_nro)) + geom_histogram(binwidth=20000)
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





