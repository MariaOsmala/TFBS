library("GenomicRanges")
library("readr")
library("data.table")
library(vroom)
library(parallel)

library(AnnotationHub)
library(GenomicFeatures)


hub <- AnnotationHub()
df=mcols(hub)
df[1:5,]
tmp=df[(df$species=="Homo sapiens" & df$genome=="GRCh38" & df$sourcetype=="GTF"),]

tmp=as.data.frame(tmp)

#tmp2=tmp["AH28674",]

tmp2=tmp["AH68821",]

unique(hub$rdataclass)

test=query(hub, c("Homo sapiens", "GRCh38", "Ensembl"))

gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf
seqinfo(gtf)["1"]

seqlengths(gtf)

setwd("/scratch/project_2006203/TFBS/")

dir()


MOODS_path="Results/MOODS2/MOODS/"
MOODS_results=dir(MOODS_path)

MOODS_file=MOODS_results[1]


# detectCores(all.tests = FALSE, logical = TRUE)

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

library("GenomicRanges")

PWMs=unique(tb$PWM)

#MOODS output is 0-based, GRanges is 1-based

motif_matches<-GRangesList()
motif_matches_smaller<-GRangesList()

for(pwm in PWMs){
  #pwm=PWMs[1]
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
  
  #choose 800,000 with the best score
  
  tmp_GRanges=tmp_GRanges[order(tmp_GRanges$score, decreasing=TRUE),]
  
  if(length(tmp_GRanges)>300000){
    tmp_GRanges_smaller=tmp_GRanges[1:300000]
    
  }
  
  motif_matches[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges
  motif_matches_smaller[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges_smaller
  
  export.bed(tmp_GRanges_smaller, paste0("Results/MOODS_bigbed/", strsplit(pwm, ".pfm")[[1]], ".bed"))
  
    
}

saveRDS(motif_matches, file = paste0("Results/MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB

saveRDS(motif_matches_smaller, file = paste0("Results/MOODS_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], "_smaller.Rds") ) #30MB







#MOODS_results=dir("/scratch/project_2006472/MOODS/")

#system("gunzip /scratch/project_2006472/MOODS/MOODS_1.csv.gz")

#Data = read.table(gzfile("MCT.txt.gz"),sep="\t")

#table=read.table(gzfile( paste0("../../Results/MOODS/", MOODS_results[1])),sep="," )

table=read.table(gzfile( paste0("/scratch/project_2006472/MOODS/", MOODS_results[1])),sep="," )

table=read.table("/scratch/project_2006472/MOODS/MOODS_1.csv",sep="," )

allData <- fread("/scratch/project_2006472/MOODS/MOODS_1.csv", showProgress = TRUE)

allData <- fread(paste0("/scratch/project_2006472/MOODS/", MOODS_results[1]), showProgress = TRUE)

zz=gzfile('file.csv.gz','rt')  
dat=read.csv(zz,header=F)
?gzfile

library(vroom)
tb <- vroom(file = paste0("../../Results/MOODS/", MOODS_results[1]), col_names = FALSE, num_threads = 25, altrep = TRUE)
tb <- vroom(file = paste0("/scratch/project_2006472/MOODS/MOODS_1.csv"), col_names = FALSE, num_threads = 25, altrep = TRUE)

tb=tb[,-7]

names(tb)=c("chr", "PWM", "position", "strand", "score", "sequence")

PWMs=unique(tb$PWM)

for(pwm in PWMs){
  
  tmp=tb[tb$PWM==pwm,]
  GRanges
  
}

#20_854_377_447
#21 GB

library(readr)
library(arrow)

#setwd("/home/local/osmalama/projects/TFBS/Results/MOODS")

data_path="/scratch/project_2006203/TFBS/Results/MOODS/"
data_path="/scratch/project_2006472/MOODS/"
fyl <- paste0(data_path,"MOODS_0.csv.gz")
pqFolder <- "/scratch/project_2006203/TFBS/Results/MOODS_compressed/"

tmp=read_csv_chunked(fyl, chunk_size=300000)

tmp=read.table(fyl,skip=0 ,nrows=300000 )



f <- function(x, pos){
  write_parquet(x,
                file.path(pqFolder, paste0(pos, ".parquet")),
                compression = "gzip",
                compression_level = 9)
}
#c = character
#i = integer
#n = number
#d = double
#l = logical

#chr1,ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm,11864,+,4.64187166754,TTTCGTTAACTTG,
read_csv_chunked(
  fyl,
  col_names=FALSE,
  callback = SideEffectChunkCallback$new(f),
  chunk_size = 10^6)


#Databases

library(inborutils)
install.packages("inborutils")

sqlite_file <- "MOODS.sqlite"
table_name <- "TF"

csv_file=paste0("../../Results/MOODS/", MOODS_results[1])

inborutils::csv_to_sqlite(csv_file = csv.name, 
                          sqlite_file, table_name, pre_process_size = 1000, 
                          chunk_size = 50000, show_progress_bar = TRUE)


#the SQlite database is available to query, similar to the previous examples:
  
  my_db <- src_sqlite("MOODS.sqlite", create = FALSE, col_names=c("chr", "name", "position", "strand", "score", "sequence"))

  TFBSs <- tbl(my_db, "TF")
  results <- TF %>%
    filter(device_info_serial == 860) %>%
    select(date_time, latitude, longitude, altitude) %>%
    filter(date_time < "2014-07-01") %>%
    filter(date_time > "2014-03-01")
  head(results)
