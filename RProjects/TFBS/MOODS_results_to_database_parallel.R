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
#devtools::install_github("cboettig/arkdb")
library("arkdb")
library("future.apply")

plan(multisession)

#Import a list of compressed tabular files (i.e. *.csv.bz2) into a local SQLite database:
  
files <- fs::dir_ls("/scratch/project_2006203/TFBS/Results/MOODS_Vierstra/")
new_db <- DBI::dbConnect(RSQLite::SQLite(), fs::path("Results/MOODS_Vierstra_SQLite/motif-hits-parallel.sqlite"))


# Strategy 2: Parallel over chunks of a table

unark(files,
  new_db, 
  streamable_table = streamable_parquet(), # required for window-parallel
  lines = 50000, 
  method = "window-parallel"
)

#Using callbacks
#It is possible to use a callback to perform just-in-time data t
#ransformations before ark writes your data object to disk in your preferred format.

unark(files, new_db, lines = 50000)


unark(files, new_db, lines = 50000)


arrays <- as.numeric(commandArgs(trailingOnly = TRUE))
#arrays=0-39

setwd("/scratch/project_2006203/TFBS/")

#results_path="/scratch/project_2006472/MOODS/"
results_path="/scratch/project_2006203/TFBS/Results/MOODS_Vierstra_processed/"

print(detectCores(all.tests = FALSE, logical = TRUE))

#hub <- AnnotationHub()
#gtf <- hub[["AH28674"]] #Homo_sapiens.GRCh38.95.gtf

#saveRDS(gtf,"RProjects/TFBS/gtf.Rds")


gtf<-readRDS("RProjects/TFBS/gtf.Rds")

con <- DBI::dbConnect(RSQLite::SQLite(), "Results/MOODS_Vierstra_SQLite/motif-hits.sqlite")
MOODS_path="/scratch/project_2006203/TFBS/Results/MOODS_Vierstra/"

for(arrays in 1:39){

start_ind=arrays*10 #100
end_ind=(arrays+1)*10-1 #100
len=10 #100

if(end_ind>398){
  end_ind=398
}


#DBI::dbDisconnect(con)


#MOODS_results=dir(MOODS_path)

for(index in seq(start_ind, end_ind, 1)){

  print(paste0("index: ", index))
  
  MOODS_file=paste0("MOODS_", index, ".csv.gz")
  
  
  
  system(paste0("gunzip ",MOODS_path,MOODS_file))
  
  
  tb <- vroom(file = paste0(MOODS_path, strsplit(MOODS_file, ".gz")[[1]]),
              col_names = FALSE, num_threads = 40, altrep = TRUE)
  
  system(paste0("gzip ",MOODS_path, strsplit(MOODS_file, ".gz")[[1]] ))
  
  
  tb=tb[,-7]
  
  names(tb)=c("chr", "motif_id", "position", "strand", "score", "sequence")
  
  
  PWMs=unique(tb$motif_id)
  
  motifs=data.frame(PWMs=do.call(rbind, strsplit(PWMs, ".pfm")), motif_id=seq(arrays*100+index*10, arrays*100+index*10 + 10-1,1))
  
  level_key=seq(arrays*100+index*10, arrays*100+index*10 + 10-1,1)
  names(level_key)=PWMs
  
  tb=tb %>% mutate(motif_id=recode(motif_id,!!!level_key))
  
  if(DBI::dbExistsTable(con, "motif-hits")){
    #TRUE
    DBI::dbAppendTable(con, "motif-hits", tb)
  }else{
    #FALSE
    DBI::dbWriteTable(con, "motif-hits", tb )
  }
  
  if(DBI::dbExistsTable(con, "motifs")){
    #TRUE
    DBI::dbAppendTable(con, "motifs", motifs)
  }else{
    #FALSE
    DBI::dbWriteTable(con, "motifs", motifs )
  }

  
  
  #MOODS output is 0-based, GRanges is 1-based
  
 # motif_matches<-GRangesList()
  
  
  
  
  
  # for(pwm in PWMs){
  #   #pwm=PWMs[1]
  #   print(pwm)
  #   tmp=tb[tb$PWM==pwm,]
  #   tmp_GRanges=GRanges(seqnames = Rle( tmp$chr, rep(1, nrow(tmp)) ),
  #   ranges = IRanges(start=tmp$position+1, end = tmp$position+ nchar(tmp$sequence[1]) ),
  #   strand = Rle(strand(tmp$strand  ), rep(1, nrow(tmp))  ))
  #   rm_index=which(names(tmp) %in% c("chr", "PWM", "position", "strand"))
  #   elementMetadata(tmp_GRanges)=tmp[, -rm_index]
  # 
  #   isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  # 
  #   seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
  #   names(seqlengths)=paste0("chr", names(seqlengths))
  # 
  #   seqlengths(tmp_GRanges)=seqlengths
  #   genome(tmp_GRanges) <- "GRCh38"
  # 
  #   #choose 300,000 with the best score
  # 
  #   tmp_GRanges=tmp_GRanges[order(tmp_GRanges$score, decreasing=TRUE),]
  # 
  #   motif_matches[[ strsplit(pwm, ".pfm")[[1]] ]]=tmp_GRanges
  #   
  # 
  #   export.bed(tmp_GRanges, paste0(results_path, "MOODS_Vierstra_bigbed/", strsplit(pwm, ".pfm")[[1]], ".bed"))
  # 
  # 
  # }
  # 
  # saveRDS(motif_matches, file = paste0(results_path, "MOODS_Vierstra_RDS/", strsplit(MOODS_file, ".csv.gz")[[1]], ".Rds")) #1.5GB
  
}
}

DBI::dbDisconnect(con)

#One could also try to store these as parquet files
