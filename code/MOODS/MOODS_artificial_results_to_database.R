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

arrays <- as.numeric(commandArgs(trailingOnly = TRUE)) #0-786

#arrays=arrays+1000

print(arrays)
setwd("/scratch/project_2006203/TFBS/")

representatives=read_tsv(file="/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives.tsv",col_names=TRUE) #3933

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

MOODS_file=dir("/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_version2.2_correct") #784 files

MOODS_file=paste0("MOODS_",arrays, ".csv.gz")

#results_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_processed/"
#results_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_union_processed/"
results_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_version2.2_processed/"





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


#MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final/"
MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_version2.2_correct/"

path_local_scratch=Sys.getenv("LOCAL_SCRATCH")
system(paste0("cp ",MOODS_path,MOODS_file, " ",path_local_scratch, "/"))

dir(path_local_scratch)
  
#system(paste0("gunzip ",path_local_scratch,"/",MOODS_file))
  
tb=read_csv(file = gzfile(paste0(path_local_scratch, "/",MOODS_file)),
            col_names = FALSE) 

#tb=read_csv(file = gzfile(paste0(MOODS_path, "/",MOODS_file)), col_names = FALSE) 
 
#tb <- vroom(file = paste0(path_local_scratch, strsplit(MOODS_file, ".gz")[[1]]),
#              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
#system(paste0("rm ",path_local_scratch, strsplit(MOODS_file, ".gz")[[1]] ))
  
  
tb=tb[,-7]
  
names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")

#There are artificial motifs for different true motifs, figure out the true motifs and the corresponding artificial  

#Need to have a list where names are the true motifs and list elements are the artificial motifs

#artificial_version2.2 <- readRDS("/projappl/project_2006203/TFBS/RProjects/TFBS/RData/artificial_version2.2.Rds")

#artificial_names <- lapply(artificial_version2.2, names)
#saveRDS(artificial_names, "/projappl/project_2006203/TFBS/RProjects/TFBS/RData/artificial_names.Rds")

artificial_names <- readRDS("/projappl/project_2006203/TFBS/RProjects/TFBS/RData/artificial_names.Rds")

PWMs_all=unique(tb$motif_id) #111
PWMs_all=gsub(".pfm", "", PWMs_all)

true_motifs=names(artificial_names)[apply(sapply(PWMs_all, 
                                                    function(y) sapply(artificial_names, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ] #111


unique_true_motifs=unique(true_motifs) #5





#motifs=data.frame(PWMs=do.call(rbind, strsplit(PWMs, ".pfm")) )

motif_matches_part=list()

for(true_motif in unique_true_motifs){
  #true_motif=unique_true_motifs[1]
  motif_matches_top<-GRangesList() 
  print(true_motif)
  #which artificial motifs correspond to this true motifs
  PWMs=paste0(PWMs_all[which(true_motifs==true_motif)], ".pfm")
  
  for(pwm in PWMs){
    # pwm=PWMs[1]
    #print(pwm)
    #pwm="ZSCAN16_HT-SELEX_TCCACC40NTTA_KR_NANTGTTAACAGAGCCTCN_2_4_YES.pfm"
    #level_key[pwm]
    tmp=tb[tb$motif_id==pwm,] #5402
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




    motif_matches_top[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_top
    
    #Save individual artificial motif bed file
    export.bed(tmp_GRanges_top, paste0(path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed"))
  
    #print(paste0("cp ", path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed ", results_path, "MOODS_bigbed/"))
    system(paste0("cp ", path_local_scratch, "/", strsplit(pwm, ".pfm")[[1]], "_top.bed ", results_path, "MOODS_bigbed/"))
  
  } #Loop over artificial
  
  motif_matches_part[[true_motif]]=motif_matches_top
  
  #Same true motif specific GRangesList
  saveRDS(motif_matches_top, file = paste0(path_local_scratch, "/", true_motif, "_top.Rds")) #1.5GB
  #paste0("cp ", path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds ", results_path, "MOODS_RDS/")
  system(paste0("cp ", path_local_scratch, "/", true_motif, "_top.Rds ", results_path, "MOODS_RDS/"))
  
  

} #Loop over true motifs

#Save arrays specific file

saveRDS(motif_matches_part, file = paste0(path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds")) #1.5GB
system(paste0("cp ", path_local_scratch, "/", strsplit(MOODS_file, ".csv.gz")[[1]], "_top.Rds ", results_path, "MOODS_RDS/"))
