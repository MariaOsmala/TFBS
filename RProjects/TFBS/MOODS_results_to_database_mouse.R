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


#arrays <- as.numeric(commandArgs(trailingOnly = TRUE))
#arrays=0-39

setwd("/scratch/project_2006203/TFBS/")

#results_path="/scratch/project_2006472/MOODS/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Vierstra_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_4_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_3_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_2_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_Mouse_processed/"
results_path="/scratch/project_2006203/TFBS/Results/MOODS_Mouse_processed_4/"
dir.create(file.path(results_path), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_bigbed"), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_RDS"), showWarnings = FALSE)

#print(detectCores(all.tests = FALSE, logical = TRUE))

#hub <- AnnotationHub()
#gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf

#saveRDS(gtf,"RProjects/TFBS/gtf.Rds")




#Uncomment if writing to database
#con <- DBI::dbConnect(RSQLite::SQLite(), "Results/MOODS_Vierstra_SQLite/motif-hits.sqlite")

#con <- DBI::dbConnect(RSQLite::SQLite(), "Results/MOODS_Teemu_SQLite/motif-hits.sqlite")

#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Vierstra/"
#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu/"
#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_4/"
#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_3/"
#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_2/"

#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39/"
MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_4/"

chrom_lengths=read.table("/projappl/project_2006203/Genomes/Mus_musculus/GRCm39_mm39/mm39.chrom.sizes")
seqlengths=chrom_lengths[,2]
names(seqlengths)=chrom_lengths[,1]
seqlengths=seqlengths[paste0("chr", c(as.character(seq(1,19,1)), "X", "Y"))]


for(arrays in 0:10){ #0:39 0:10 0:7 0:5

start_ind=arrays*10 #100
end_ind=(arrays+1)*10-1 #100
len=10 #100

#3982+15=3997

#Old
#if(end_ind>398){ # 398 
#  end_ind=398
#}

#if(end_ind>399){ # 398 
#  end_ind=399
#}

if(end_ind>100){ #  102
  end_ind=100
}

#if(end_ind>76){ #  76
#  end_ind=76
#}

#if(end_ind>57){ #  57
#  end_ind=57
#}


#DBI::dbDisconnect(con)


#MOODS_results=dir(MOODS_path)

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
  

  
  
  motif_matches<-GRangesList()
   motif_matches_top<-GRangesList()
   
   for(pwm in PWMs){
   #   #pwm=PWMs[1]
      print(pwm)
     #pwm="ZSCAN16_HT-SELEX_TCCACC40NTTA_KR_NANTGTTAACAGAGCCTCN_2_4_YES.pfm"
      #level_key[pwm]
      tmp=tb[tb$motif_id==level_key[pwm],]
      tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
      ranges = IRanges(start=tmp$position+1, end = tmp$position+ nchar(tmp$sequence[1]) ),
      strand = Rle(strand(tmp$strand  ), rep(1, nrow(tmp))  ))
      rm_index=which(names(tmp) %in% c("chr", "motif_id", "position", "strand"))
      elementMetadata(tmp_GRanges)=tmp[, -rm_index]
   #
      isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
   #
      
      
   #
  
      #names(seqlengths) %in% names(seqlengths(tmp_GRanges))
      seqlengths(tmp_GRanges)=seqlengths[names(seqlengths) %in% names(seqlengths(tmp_GRanges))]
      genome(tmp_GRanges) <- "mm39"
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
      motif_matches[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges
      motif_matches_top[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_top
   #
   #
      export.bed(tmp_GRanges, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], ".bed"))
      export.bed(tmp_GRanges_top, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], "_top.bed"))
   #
   #
    }
   #
    saveRDS(motif_matches, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB
    saveRDS(motif_matches_top, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds")) #1.5GB

}
} #0-39


#DBI::dbDisconnect(con)

#One could also try to store these as parquet files
