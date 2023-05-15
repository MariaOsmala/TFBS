library("GenomicRanges")
library("readr")
library("data.table")
library(vroom)
library(parallel)
library(AnnotationHub)
library(GenomicFeatures)
library(rtracklayer)


arrays <- as.numeric(commandArgs(trailingOnly = TRUE))

setwd("/scratch/project_2006203/TFBS/")

results_path="/scratch/project_2006472/MOODS/"


print(detectCores(all.tests = FALSE, logical = TRUE))

#hub <- AnnotationHub()
#gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf

#saveRDS(gtf,"RProjects/TFBS/gtf.Rds")

gtf<-readRDS("RProjects/TFBS/gtf.Rds")
#/scratch/project_2006203/TFBS/RProjects/TFBS/gtf.Rds

start_ind=arrays*10 #100
end_ind=(arrays+1)*10-1 #100
len=10 #100

if(end_ind>398){
  end_ind=398
}


MOODS_path="Results/MOODS2/MOODS/"
#MOODS_results=dir(MOODS_path)

for(index in seq(start_ind, end_ind, 1)){

  print(paste0("index: ", index))
  
  MOODS_file=paste0("MOODS_", index, ".csv.gz")
  
  
  
  system(paste0("gunzip ",MOODS_path,MOODS_file))
  
  
  tb <- vroom(file = paste0(MOODS_path, strsplit(MOODS_file, ".gz")[[1]]),
              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
  system(paste0("gzip ",MOODS_path, strsplit(MOODS_file, ".gz")[[1]] ))
  
  
  tb=tb[,-7]
  
  names(tb)=c("chr", "PWM", "position", "strand", "score", "sequence")
  
  #chr: Sequence name (either file name or fasta header).
  #PWM: Matrix file name.
  #position: Hit position (first position of sequence is position 0).
  #strand: Indicates whether the match is versus the input strand (+) or the reverse complement (-).
  #score: Match score.
  #sequence: Hit site in the original sequence (note that if the match is -, then the match is versus the reverse complement of this).
  #7th column: Hit sequence with SNPs applied, i.e. shows how the ambiguity codes are resolved to get this score
  #(there can be multiple hits at the same position with different "real" hit sequences).
  
  
  #chr   PWM                                                            position strand score sequence
  #<chr> <chr>                                                             <dbl> <chr>  <dbl> <chr>
  #1 chr1  ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm    11864 +       4.64 TTTCGTTAACTTG
  #2 chr1  ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm    12144 -       5.01 TAACTTAATACCA
  #3 chr1  ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm    12563 +       4.33 AGCTGGTGATGTG
  #4 chr1  ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm    16361 +       5.04 TTGGTTTAATTAG
  #5 chr1  ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm    16509 +       4.42 AATTGAGAATTTC
  
  #nchar(tmp$sequence[1])
  
  #tmp$position[2]+1
  #tmp$position[2]+nchar(tmp$sequence[2])
  
  #TAACTTAATACCA #sequence 2
  #ATTGAATTATGGT # complement
  
  #TTTCGTTAACTTG #first seqquence
  #TGGTATTAAGTTA #reverse complement
  #AGCTGGTGATGTG #third
  #TTGGTTTAATTAG #fourth
  #AATTGAGAATTTC #fifth
  
  
  
  PWMs=unique(tb$PWM)
  
  #MOODS output is 0-based, GRanges is 1-based
  
  motif_matches<-GRangesList()
  motif_matches_smaller<-GRangesList()
  
  for(pwm in PWMs){
    #pwm=PWMs[1]
    print(pwm)
    tmp=tb[tb$PWM==pwm,]
    tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
    ranges = IRanges(start=tmp$position+1, end = tmp$position+ nchar(tmp$sequence[1]) ),
    strand = Rle(strand(tmp$strand  ), rep(1, nrow(tmp))  ))
    rm_index=which(names(tmp) %in% c("chr", "PWM", "position", "strand"))
    elementMetadata(tmp_GRanges)=tmp[, -rm_index]
  
    isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  
    seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
    names(seqlengths)=paste0("chr", names(seqlengths))
  
    seqlengths(tmp_GRanges)=seqlengths
    genome(tmp_GRanges) <- "GRCh38"
  
    #choose 300,000 with the best score
  
    tmp_GRanges=tmp_GRanges[order(tmp_GRanges$score, decreasing=TRUE),]
  
    if(length(tmp_GRanges)>300000){
      tmp_GRanges_smaller=tmp_GRanges[1:300000]
  
    }
  
    motif_matches[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges
    motif_matches_smaller[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_smaller
  
    export.bed(tmp_GRanges_smaller, paste0(results_path, "MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], ".bed"))
  
  
  }
  
  saveRDS(motif_matches, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB
  
  saveRDS(motif_matches_smaller, file = paste0(results_path, "MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], "_smaller.Rds") ) #30MB


}
#One could also try to store these as parquet files
