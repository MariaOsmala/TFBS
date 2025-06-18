#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#or 
#lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))




library("GenomicRanges")
library(vroom)
library(rtracklayer)
library("GenomeInfoDb")

print(sessionInfo())
print(.libPaths())

args <- commandArgs(trailingOnly = TRUE)
arrays=as.numeric(args[1])
MOODS_path=args[2]
#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset/"
results_path=args[3]
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset_processed/" 

#setwd("/scratch/project_2006203/TFBS/")

dir.create(file.path(results_path), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_bigbed"), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_RDS"), showWarnings = FALSE)

seqinfo_obj <- Seqinfo(genome="hg38")  # or appropriate genome build

#for(arrays in 0:33){ #0:39 0:10 0:7 0:5

start_ind=arrays*10 #100
end_ind=(arrays+1)*10-1 #100
len=10 #100

 #3933 motifs
if(end_ind>393){ # 398 
 end_ind=393
}
 
 

for(index in seq(start_ind, end_ind, 1)){ #0-9
  #index=0
  print(paste0("index: ", index))
  
  MOODS_file=paste0("MOODS_", index, ".csv.gz")
  
  
  
  system(paste0("gunzip ",MOODS_path,MOODS_file))
  
  
  tb <- vroom(file = paste0(MOODS_path, strsplit(MOODS_file, ".gz")[[1]]),
              col_names = FALSE)
  
  system(paste0("gzip ",MOODS_path, strsplit(MOODS_file, ".gz")[[1]] ))
  
  
  tb=tb[,-7]
  
  names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")
  
  
  PWMs=unique(tb$motif_id)

  
  #MOODS output is 0-based, GRanges is 1-based
  
  #motif_matches<-GRangesList()
  motif_matches_top<-GRangesList()
   
   for(pwm in PWMs){
   #   #pwm=PWMs[1]
      print(pwm)
      tmp=tb[tb$motif_id==pwm,]
      tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
      ranges = IRanges(start=tmp$position+1, end = tmp$position+nchar(tmp$sequence[1]) ),
      strand = Rle(strand(tmp$strand  ), rep(1, nrow(tmp))  ))
      rm_index=which(names(tmp) %in% c("chr", "motif_id", "position", "strand"))
      elementMetadata(tmp_GRanges)=tmp[, -rm_index]

      
      seqinfo(tmp_GRanges)=seqinfo_obj[seqlevels(tmp_GRanges)]
      
      genome(tmp_GRanges) <- "hg38"
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
  
  
  
   #
      #motif_matches[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges
      motif_matches_top[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_top
   #
   #
      #export.bed(tmp_GRanges, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], ".bed"))
      export.bed(tmp_GRanges_top, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], "_top.bed"))
   #
   #
    }
   #
    #saveRDS(motif_matches, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB
    saveRDS(motif_matches_top, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds")) #1.5GB

}
#} #0-39



