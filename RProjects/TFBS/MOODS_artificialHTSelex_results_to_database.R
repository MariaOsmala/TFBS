#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#or 
#lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))


print(sessionInfo())
print(.libPaths())

library("GenomicRanges")
library("readr")
library("data.table")
library(vroom)
library(parallel)
library(AnnotationHub)
library(GenomicFeatures)
library(rtracklayer)
library("dbplyr")
library("dplyr")
library(BSgenome.Hsapiens.UCSC.hg38)

print(sessionInfo())
print(.libPaths())


setwd("/scratch/project_2006203/TFBS/")



results_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final_artificialHTSelex_processed/"
dir.create(file.path(results_path), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_bigbed"), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_RDS"), showWarnings = FALSE)


gtf<-readRDS("RProjects/TFBS/gtf.Rds")


MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final_artificialHTSelex/" 

MOODS_file=paste0("MOODS.csv.gz")
#system(paste0("gunzip ",MOODS_path,MOODS_file))
  
  
tb <- vroom(file = paste0(MOODS_path, strsplit(MOODS_file, ".gz")[[1]]),
              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
#system(paste0("gzip ",MOODS_path, strsplit(MOODS_file, ".gz")[[1]] ))
  
tb=tb[,-7]
names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")
  
PWMs=unique(tb$motif_id)
motifs=data.frame(PWMs=do.call(rbind, strsplit(PWMs, ".pfm")) )
  
#MOODS output is 0-based, GRanges is 1-based
  
#motif_matches<-GRangesList()
motif_matches_top<-GRangesList()
   
 for(pwm in PWMs){
   #pwm=PWMs[1]
    print(pwm)
    tmp=tb[tb$motif_id==pwm,]
    tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
    ranges = IRanges(start=tmp$position+1, end = tmp$position+nchar(tmp$sequence[1]) ),
    strand = Rle(strand(tmp$strand  ), rep(1, nrow(tmp))  ))
    rm_index=which(names(tmp) %in% c("chr", "motif_id", "position", "strand"))
    elementMetadata(tmp_GRanges)=tmp[, -rm_index]
 #
    isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
 #
    seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
    names(seqlengths)=paste0("chr", names(seqlengths))
 #

    #names(seqlengths) %in% names(seqlengths(tmp_GRanges))
    seqlengths(tmp_GRanges)=seqlengths[names(seqlengths) %in% names(seqlengths(tmp_GRanges))]
    genome(tmp_GRanges) <- "GRCh38"
 #
 #   #choose 300,000 with the best score
 #
    tmp_GRanges=tmp_GRanges[order(tmp_GRanges$score, decreasing=TRUE),]
#
    #what is the smallest score in the top 300,000
    if(length(tmp_GRanges)>300000){
     score_min=tmp_GRanges[300000]$score
     tmp_GRanges_top=tmp_GRanges[tmp_GRanges$score >= score_min]
    }else{
      tmp_GRanges_top=tmp_GRanges
    }



    motif_matches_top[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_top

    export.bed(tmp_GRanges_top, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], "_top.bed"))

  }


saveRDS(motif_matches_top, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds")) #1.5GB





