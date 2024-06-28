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

arrays <- as.numeric(commandArgs(trailingOnly = TRUE))
#arrays=0-39
print(arrays)
setwd("/scratch/project_2006203/TFBS/")


#results_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final_union_processed/"
results_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct_processed/"
dir.create(file.path(results_path), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_bigbed"), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_RDS"), showWarnings = FALSE)

#print(detectCores(all.tests = FALSE, logical = TRUE))

#hub <- AnnotationHub()
#gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf

#saveRDS(gtf,"RProjects/TFBS/gtf.Rds")

#Is this the right genome version, YES, the genome lengths of gtf and bsg from Bsgenome are the same
gtf<-readRDS("RProjects/TFBS/gtf.Rds")

#1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16 
#248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189  90338345 
#17        18        19        20        21        22         X         Y 
#83257441  80373285  58617616  64444167  46709983  50818468 156040895  57227415 

#
#bsg <- BSgenome.Hsapiens.UCSC.hg38
#seq_len_GRCh38=seqlengths(bsg)[paste0("chr",c(as.character(seq(1,22,1)), "X", "Y"))]
#chr1      chr2      chr3      chr4      chr5      chr6      chr7      chr8      chr9     chr10     chr11     chr12     chr13     chr14     chr15     chr16 
#248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189  90338345 
#chr17     chr18     chr19     chr20     chr21     chr22      chrX      chrY 
#83257441  80373285  58617616  64444167  46709983  50818468 156040895  57227415 


#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_human_final/"
#MOODS_path="/scratch/project_2006203/TFBS/Results//MOODS_human_final_union/" 
MOODS_path="/scratch/project_2006203/TFBS/Results//MOODS_human_final_version2.2_correct/" 

#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39/"


#for(arrays in 0:33){ #0:39 0:10 0:7 0:5

 start_ind=arrays*10 #100
 end_ind=(arrays+1)*10-1 #100
 len=10 #100

#start_ind=arrays
#end_ind=arrays

#3982+15=3997

#Old
#if(end_ind>398){ # 398 
#  end_ind=398
#}

#if(end_ind>399){ # 398 
#  end_ind=399
#}

#3294 motifs

#if(end_ind>329){ # 398 
#  end_ind=329
#}

#4003 motifs
#if(end_ind>400){ # 398 
#  end_ind=400
#}

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
              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
  system(paste0("gzip ",MOODS_path, strsplit(MOODS_file, ".gz")[[1]] ))
  
  
  tb=tb[,-7]
  
  names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")
  
  
  PWMs=unique(tb$motif_id)
  
  #motif_id=seq(arrays*100+index, arrays*100+index + 10-1,1)
  
  #motif_id=seq(index*10, index*10 + 3-1,1)
  
  #motif id starts from zero
  
  #Uncomment if using database ######################
  motifs=data.frame(PWMs=do.call(rbind, strsplit(PWMs, ".pfm")), motif_id=seq(index*10, index*10 + length(PWMs)-1,1) )
  level_key=seq(index*10, index*10 + length(PWMs)-1,1)
  names(level_key)=PWMs


  tb=tb %>% mutate(motif_id=recode(motif_id,!!!level_key))
  
  #MOODS output is 0-based, GRanges is 1-based
  
  #motif_matches<-GRangesList()
  motif_matches_top<-GRangesList()
   
   for(pwm in PWMs){
   #   #pwm=PWMs[1]
      print(pwm)
      tmp=tb[tb$motif_id==level_key[pwm],]
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



