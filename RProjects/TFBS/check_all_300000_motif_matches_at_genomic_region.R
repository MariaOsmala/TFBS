library("rtracklayer")
library("TFBSTools")
library("Biostrings")
library("motifmatchr")
library("GenomicRanges")

#.libPaths(c("/projappl/project_2006203/project_rpackages", .libPaths()))


#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames_match_numbers_MOODS_threshold.tsv", 
                           sep="\t", header=TRUE)

rep_motifs=representatives[which(representatives$new_representative=="YES"),]


#https://www.ncbi.nlm.nih.gov/gene/105804841

#chr7:156,790,708-156,793,079 #this is 2,372 bp.

#coordinates in GRCh37 chr7:156583402..156585773


#Vista browser, this is in hg19:  chr7:156,583,782-156,584,569


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
genome(test_GRanges_GRCh38) <- "GRCh38"

Vista_GRanges_GRCh38<- test_GRanges_GRCh38

export.bed(test_GRanges_GRCh38, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider/VISTA_SHH_enhancer_hs2496.bed"))



#Wider region 2372

#This is hg19
test_GRanges=GRanges(seqnames = Rle( "chr7", rep(1, 1) ),
                     ranges = IRanges(start=156583402, end = 156585773) )


#GRCh37/hg19 

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



for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  #i=1
  threshold=rep_motifs$MOODS_threshold[i]
  name=rep_motifs$ID[i]
  #rep_motifs$match_numbers[i]
  if(threshold==5){
    TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/",name,"_top.bed"))
  }else{
    TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_",threshold,"_processed/MOODS_bigbed/",name,"_top.bed"))
  }
  
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
  export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider/bed/",rep_motifs$ID[i],".bed") )
  }

}

#Collect to a table

GRanges_all<-GRanges()

for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  #i=1
  threshold=rep_motifs$MOODS_threshold[i]
  name=rep_motifs$ID[i]
  #rep_motifs$match_numbers[i]
  if(threshold==5){
    TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/",name,"_top.bed"))
  }else{
    TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_",threshold,"_processed/MOODS_bigbed/",name,"_top.bed"))
  }
  
  isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  
  seqlengths=seqlengths(seqinfo(gtf_GRCh38))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  
  seqlevels(TFBS_GRanges)=paste0("chr",seqlevels(seqinfo(keepSeqlevels(gtf_GRCh38, c(as.character(seq(1,22,1)), "X", "Y"), pruning.mode="coarse"))))
  
  seqlengths(TFBS_GRanges)=seqlengths
  
  genome(TFBS_GRanges) <- "GRCh38"
  
  print(name)
  print(length(findOverlaps( TFBS_GRanges,Vista_GRanges_GRCh38, ignore.strand=TRUE)))
  
  hits=findOverlaps( TFBS_GRanges,Vista_GRanges_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
    
    tmp=TFBS_GRanges[hits@from]
    tmp$name=name
    GRanges_all=c(GRanges_all, tmp)
  }
  
}


#Let's sort GRanges_all by score

tmp=as.data.frame(GRanges_all[order(GRanges_all$score, decreasing=TRUE)])
write.table(tmp, quote=FALSE,file="/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider/SHH_enhancer_ZRS_matches.tsv",sep="\t", row.names = FALSE)
