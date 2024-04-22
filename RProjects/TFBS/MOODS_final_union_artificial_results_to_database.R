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

print(sessionInfo())
print(.libPaths())
.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")

arrays <- as.numeric(commandArgs(trailingOnly = TRUE)) #644

#arrays=arrays+1000

print(arrays)
setwd("/scratch/project_2006203/TFBS/")

representatives=read_csv(file="/projappl/project_2006203/TFBS/PWMs_final/representatives.csv",col_names=FALSE)

#MOODS_PROX1_HOXA2_TCATAA40NAGC_YJI_NTAATTAAAGAN_m1_c3b0_short_composite_new.csv.gz
#which(representatives$X1=="ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420") 767

#df.representatives=read.table("/projappl/project_2006203/TFBS/PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
#df.representatives=df.representatives[which(df.representatives$new_representative=="YES"),]


#MOODS_BCL6B_HT-SELEX_TGCGGG20NGA_AC_TGCTTTCTAGGAATTMM_2_4.zst Does not work
#library("fst")
# read compressed file into a raw vector

#read_csv(file = gzfile(paste0(path_local_scratch, "/",MOODS_file)),
#         col_names = FALSE) 

#compressed_file="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_zstd/MOODS_BCL6B_HT-SELEX_TGCGGG20NGA_AC_TGCTTTCTAGGAATTMM_2_4.zst"

#compressed_vec <- readBin(compressed_file, "raw", file.size(compressed_file))
# decompress file contents
#raw_vec_decompressed <- decompress_fst(compressed_vec)

#A nice feature of data.table’s fread method is that 
#it can parse in-memory data directly. That means that we can easily feed our raw vector to fread:
  
#library(data.table)
# read data set from the in-memory csv
#dt <- fread(rawToChar(raw_vec_decompressed))



MOODS_file=paste0("MOODS_",representatives[arrays+1,]$X1, ".csv.gz")

results_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_processed/"



#which runs fail
#finished=dir(paste0(results_path, "MOODS_RDS/"))
#finished=gsub("MOODS_", "", finished)
#finished=gsub("_top.Rds", "", finished)

#failed=which(!(representatives$X1 %in% finished))[which(!(representatives$X1 %in% finished))<800]
#finished_id=which((representatives$X1 %in% finished))[which((representatives$X1 %in% finished))<800]

#those which fail have a low information content
#hist(df.representatives$IC[failed])
#hist(df.representatives$IC[finished_id])

dir.create(file.path(results_path), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_bigbed"), showWarnings = FALSE)
dir.create(file.path(results_path, "MOODS_RDS"), showWarnings = FALSE)

#print(detectCores(all.tests = FALSE, logical = TRUE))

#hub <- AnnotationHub()
#gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf

#saveRDS(gtf,"RProjects/TFBS/gtf.Rds")

#Is this the right genome version, YES, the genome lengths of gtf and bsg from Bsgenome are the same
gtf<-readRDS("RProjects/TFBS/gtf.Rds")

MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final/"

path_local_scratch=Sys.getenv("LOCAL_SCRATCH")
system(paste0("cp ",MOODS_path,MOODS_file, " ",path_local_scratch, "/"))

dir(path_local_scratch)
  
#system(paste0("gunzip ",path_local_scratch,"/",MOODS_file))
  
tb=read_csv(file = gzfile(paste0(path_local_scratch, "/",MOODS_file)),
            col_names = FALSE) 
 
#tb <- vroom(file = paste0(path_local_scratch, strsplit(MOODS_file, ".gz")[[1]]),
#              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
#system(paste0("rm ",path_local_scratch, strsplit(MOODS_file, ".gz")[[1]] ))

library(processx)

read_zstd_csv <- function(file_path, ...) {
  # Create a process to decompress the file using zstd
  p <- process$new("zstd", args = c("-d", "--stdout", file_path), stdout = "|")
  
  # Use the connection to read the decompressed data into R
  data <- read.csv(p$get_output_connection(), ...)
  
  # Ensure the external process is terminated
  p$kill()
  
  return(data)
}

# Example usage:
# Replace 'path/to/your/data.zst' with the actual path to your zstd compressed file
data_frame <- read_fst(path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_union_zstd/MOODS_0.zst")

data_file <- file.path("/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_union_zstd/MOODS_0.zst")
data <- readBin(data_file, raw(), file.info(data_file)$size)
decomp=zstdDecompress(data,file.info(data_file)$size)

library(fst)
  
tb=tb[,-7]
  
names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")
  
  
PWMs=unique(tb$motif_id)
  
motifs=data.frame(PWMs=do.call(rbind, strsplit(PWMs, ".pfm")) )

motif_matches_top<-GRangesList()
   
for(pwm in PWMs){
#   #pwm=PWMs[1]
  print(pwm)
 #pwm="ZSCAN16_HT-SELEX_TCCACC40NTTA_KR_NANTGTTAACAGAGCCTCN_2_4_YES.pfm"
  #level_key[pwm]
  tmp=tb[tb$motif_id==pwm,]
  tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
  ranges = IRanges(start=tmp$position+1, end = tmp$position+ nchar(tmp$sequence[1]) ),
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
  export.bed(tmp_GRanges_top, paste0(path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed"))
  #IS THERE SOMETHING WRONG?
  print(paste0("cp ", path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed ", results_path, "MOODS_bigbed/"))
  system(paste0("cp ", path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed ", results_path, "MOODS_bigbed/"))
  
  
#
#
}
#
#saveRDS(motif_matches, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB
saveRDS(motif_matches_top, file = paste0(path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds")) #1.5GB
paste0("cp ", path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds ", results_path, "MOODS_RDS/")
system(paste0("cp ", path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds ", results_path, "MOODS_RDS/"))


