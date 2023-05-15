library(readr)
library(arrow)
setwd("/home/local/osmalama/projects/TFBS/Results/MOODS")


data_path="/home/local/osmalama/csc_scratch/TFBS/Results/MOODS/"
data_path="/home/local/osmalama/csc_CRISPR_scratch/MOODS/"

MOODS_results=dir(data_path)
library(vroom)
tb <- vroom(file = "MOODS_0.csv.gz", col_names = FALSE, num_threads = 25, altrep = TRUE)

tb <- vroom(file = paste0("/home/local/osmalama/csc_CRISPR_scratch/MOODS/", MOODS_results[1]))

tb <- vroom(file = paste0("/home/local/osmalama/csc_CRISPR_scratch/MOODS/MOODS_1.csv"))
#20_854_377_447
#21 GB


unpigz -c MOODS_10.csv.gz | wc -l
#1_393_675_622

#unpigz -c MOODS_0.csv.gz | wc -l
#1_386_410_656


#wc -l MOODS_0.csv 
#1_386_410_656 MOODS_0.csv

getwd()

df <- read_parquet("1.parquet", as_data_frame=TRUE)


?read_parquet


8000*300000

fyl <- paste0(data_path,"MOODS_28.csv.gz")
fyl <- paste0(data_path,"MOODS_0.csv")

system("zcat /home/local/osmalama/csc_scratch/TFBS/Results/MOODS/MOODS_28.csv.gz | wc -l", TRUE)

system("wc -l /home/local/osmalama/csc_scratch/TFBS/Results/MOODS/MOODS_0.csv | wc -l", TRUE)

#There can be more than 300000 hits (this is just at least), not exact
#Choose the 300000 with best scores

table=read.table(fyl,skip=0 ,nrows=900000, sep="," )
length(which(table[,2]==unique(table[,2])[1] ))

table2=read.table(fyl,skip=300000 ,nrows=300000, sep="," )
table3=read.table(fyl,skip=209367368-300000 ,nrows=300000, sep="," )

209367368-300000

unique( table2[,2])
unique( table3[,2])
pqFolder <- "/home/local/osmalama/projects/TFBS/Results/MOODS/"

f <- function(x, pos){
  write_parquet(x,
                file.path(pqFolder, paste0(pos, ".parquet")),
                compression = "gzip",
                compression_level = 9)
}
#chr1,ALX3_TBX20_TCCAAA40NTCA_YWIIII_NGGTGTTAATTRN_m1_c3b0_short.pfm,11864,+,4.64187166754,TTTCGTTAACTTG,
read_csv_chunked(
  fyl,
  col_types = list(Column1="chr", Column2="PWM", Column3="position", Column4="strand", 
                   Column5="score", Column6="Hit-site", Column7="Hit-site-SNP"), # all column specifications
  callback = SideEffectChunkCallback$new(f),
  chunk_size = 10^6)
