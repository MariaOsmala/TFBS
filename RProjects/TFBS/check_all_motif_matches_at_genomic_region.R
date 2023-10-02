invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))


library("rtracklayer")
library("TFBSTools")
library("Biostrings")
library("motifmatchr")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1", .libPaths()))

.libPaths("/projappl/project_2006203/project_rpackages_4.2.1", include.site=FALSE)

.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", force=TRUE)

BiocManager::install("motifmatchr", site_repository = "/projappl/project_2006203/project_rpackages_4.2.1", force=TRUE, update=TRUE)
install.packages("/projappl/project_2006203/project_rpackages_4.2.1/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", repos = NULL, type = "source")

#echo "R_LIBS=/projappl/project_2006203/project_rpackages_4.2.1" >> ~/.Renviron

#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames_match_numbers_MOODS_threshold.tsv", 
#                           sep="\t", header=TRUE)


representatives=read.table("/projappl/project_2006203/TFBS/PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294


rep_motifs=representatives[which(representatives$new_representative=="YES"),]

#remove _pfm_composite_new and _pfm_spacing_new
#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

#https://www.ncbi.nlm.nih.gov/gene/105804841 GRCh38: chr7:156,790,708-156,793,079 #this is 2,372 bp.

#coordinates in GRCh37 chr7:156583402..156585773

#Vista browser, this is in hg19:  chr7:156,583,782-156,584,569


  
i=1


motifs<-PFMatrixList()


for(i in 1:nrow(rep_motifs)){
  print(i)
  profileMatrix=as.matrix(read.table(file=paste0("/projappl/project_2006203/TFBS/",gsub("../../", "",rep_motifs$filename[i]))))
  dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
  
motifs[[rep_motifs$ID[i]]] <- PFMatrix(ID=rep_motifs$ID[i], name=rep_motifs$symbol[i], 
                #matrixClass="Zipper-Type", 
                strand="+", 
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                tags=list(family=rep_motifs$Lambert2018.families[i], 
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

#saveRDS(motifs, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pfms.Rds")
motifs<-readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pfms.Rds")

#motifmatchr expects input PWMs to use either natural logarithm or log 2. 
#If the input is a PFM, the TFBSTools toPWM is used for making the PWM 
#with the default psueodcounts of 0.8 and base 2 logarithm. For more control of the pseudocounts, 
#simply use the toPWM function to convert your PFMs prior to calling matchMotifs.

#MOODS
#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#  --log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#  --lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

library(TFBSTools)
example_pwms <- do.call(PWMatrixList,lapply(motifs, toPWM, 
                                            pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))

#saveRDS(example_pwms, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")
example_pwms=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")

gtf_GRCh38 <- readRDS("/scratch/project_2006203/TFBS/RProjects/TFBS/gtf.Rds")

# Get motif matches for example motifs in peaks 

#The VISTA enhancer, 787 bp region
#156584569-156583782 This is hg19, convert to GRCh38
test_GRanges=GRanges(seqnames = Rle( "chr7", rep(1, 1) ),
                     ranges = IRanges(start=156583782, end = 156584569) )

genome(test_GRanges) <- "GRCh37"

ChainFilePath='/projappl/project_2006203/liftOver/hg19ToHg38.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
test_GRanges_GRCh38=unlist(rtracklayer::liftOver(test_GRanges, ch))
isCircular(seqinfo(test_GRanges_GRCh38))=rep(FALSE, length(seqinfo(test_GRanges_GRCh38)))
seqlengths=seqlengths(seqinfo(gtf_GRCh38))[as.character(7)]
names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(test_GRanges_GRCh38)=seqlengths
genome(test_GRanges_GRCh38) <- "hg38"

export.bed(test_GRanges_GRCh38, paste0("/scratch/project_2006203/TFBS/Results/SHH_enhancer_ZRS/VISTA_SHH_enhancer_hs2496.bed"))


motif_ix <- matchMotifs(example_pwms, test_GRanges_GRCh38, genome="hg38", p.cutoff=1e-05, out="positions") #out=scores


#chr1:214150857-214151185
#sort based on the max score

motif_ix_maxscore<-list()

motif_ix_lnscore<-list()

for(i in 1:length(motif_ix)){
  print(i)
  motif_ix[[i]]=motif_ix[[i]][rev(order(motif_ix[[i]]$score))]
  motif_ix_lnscore[[i]]=motif_ix[[i]]
  motif_ix_lnscore[[i]]$score=motif_ix[[i]]$score/log2(exp(1))
  motif_ix_lnscore[[i]]=motif_ix_lnscore[[i]][motif_ix_lnscore[[i]]$score>=2]
  if(length(motif_ix_lnscore[[i]]) != 0){
      export.bed(motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/SHH_enhancer_ZRS/",names(motif_ix)[[i]],".bed") )
      motif_ix_maxscore[[names(motif_ix)[[i]] ]]=as.data.frame(motif_ix_lnscore[[i]][1,])
  }
}




best_motif_matches=do.call(rbind, motif_ix_maxscore)
best_motif_matches$pfm=names(motif_ix_maxscore)

#sort based on motif match
best_motif_matches=best_motif_matches[order(best_motif_matches$score, decreasing=TRUE),]

#Wider region 2372

#This is hg19
#test_GRanges=GRanges(seqnames = Rle( "chr7", rep(1, 1) ),
#                     ranges = IRanges(start=156583402, end = 156585773) )

#In hg19 coordinates chr7:156583402..156585773



# motif_ix_2kb <- matchMotifs(example_pwms, test_GRanges, genome = "hg19", p.cutoff=1e-05, out="positions") #out=scores
# 
# motif_ix_2kb_maxscore<-list()
# 
# motif_ix_2kb_lnscore<-list()
# 
# for(i in 1:length(motif_ix_2kb)){
#   print(i)
#   motif_ix_2kb[[i]]=motif_ix_2kb[[i]][rev(order(motif_ix_2kb[[i]]$score))]
#   motif_ix_2kb_lnscore[[i]]=motif_ix_2kb[[i]]
#   motif_ix_2kb_lnscore[[i]]$score=motif_ix_2kb[[i]]$score/log2(exp(1))
#   motif_ix_2kb_lnscore[[i]]=motif_ix_2kb_lnscore[[i]][motif_ix_2kb_lnscore[[i]]$score>=2]
#   if(length(motif_ix_2kb_lnscore[[i]]) != 0){
#     export.bed(motif_ix_2kb_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/SHH_enhancer_ZRS_wider/",names(motif_ix_2kb)[[i]],".bed") )
#     motif_ix_2kb_maxscore[[names(motif_ix_2kb)[[i]] ]]=as.data.frame(motif_ix_2kb_lnscore[[i]][1,])
#   }
# }


# best_motif_matches_2kb=do.call(rbind, motif_ix_2kb_maxscore)
# best_motif_matches_2kb$pfm=names(motif_ix_2kb_maxscore)

#sort based on motif match
# best_motif_matches_2kb=best_motif_matches_2kb[order(best_motif_matches_2kb$score, decreasing=TRUE),]


#save.image(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/SHH_enhancer_ZRS.Rds")
#load(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/SHH_enhancer_ZRS.Rds")



#GRCh37/hg19 chr. 1: 214,150,857-214,151,185

#gtf_GRCh37 <- readRDS("/projappl/project_2006203/chrom_lengths/gtf_GRCh37.Rds")
#test region wider
test_GRanges=GRanges(seqnames = Rle( "chr7", rep(1, 1) ),
                     ranges = IRanges(start=156583402, end = 156585773),
                     strand = Rle(strand("*"), rep(1, 1)) 
                     )


#seqlevelsStyle(test_GRanges)="NCBI"
#seqlevels(test_GRanges)=seqlevels(seqinfo(keepSeqlevels(gtf_GRCh37, c(7), pruning.mode="coarse")))
#seqinfo(test_GRanges)=seqinfo(keepSeqlevels(gtf_GRCh37, 7, pruning.mode="coarse"))
#seqlevelsStyle(test_GRanges) #NCBI
genome(test_GRanges) <- "GRCh37"

ChainFilePath='/projappl/project_2006203/liftOver/hg19ToHg38.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)

#chr7 156,790,708-156,793,079      * #chr7:156,790,708-156,793,079 THIS IS CORRECT
test_GRanges_GRCh38=unlist(rtracklayer::liftOver(test_GRanges, ch))


isCircular(seqinfo(test_GRanges_GRCh38))=rep(FALSE, length(seqinfo(test_GRanges_GRCh38)))

seqlengths=seqlengths(seqinfo(gtf_GRCh38))[as.character(7)]
names(seqlengths)=paste0("chr", names(seqlengths))

seqlengths(test_GRanges_GRCh38)=seqlengths
genome(test_GRanges_GRCh38) <- "GRCh38"

#GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds")) #
GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")) #

for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  #print(i)
  
  #i=1
  #threshold=rep_motifs$MOODS_threshold[i]
  name=rep_motifs$ID[i]
  #rep_motifs$match_numbers[i]
  #if(threshold==5){
  #  TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/",name,"_top.bed"))
  #}else{
  #  TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_",threshold,"_processed/MOODS_bigbed/",name,"_top.bed"))
  #}
  
  TFBS_GRanges=GR_list[[name]]
  
  isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  
  seqlengths=seqlengths(seqinfo(gtf_GRCh38))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  
  seqlevels(TFBS_GRanges)=paste0("chr",seqlevels(seqinfo(keepSeqlevels(gtf_GRCh38, c(as.character(seq(1,22,1)), "X", "Y"), pruning.mode="coarse"))))
  
  seqlengths(TFBS_GRanges)=seqlengths
  
  genome(TFBS_GRanges) <- "GRCh38"
  
  print(name)
  print(length(findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE)))
  
  hits=findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
  export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/",rep_motifs$ID[i],".bed") )
  }

}












