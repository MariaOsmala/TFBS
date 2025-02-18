# Libraries and parameters --------------------------------------------------------------------


library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(httr)
library(AnnotationHub)
library(tidyverse)

data_path="../RData"

# Extract info, coordinates, and non-coding regions of protein coding genes --------------------------------------------------------------------

#Extract the canonical transcripts, their TSS, and 5' utr regions, also exons (needed in CRISPR-screening data analysis)
#Only in the canonical chromosomes
#Query also exons

## Chromosome lengths for GRCh38 and GRCh37 --------------------------------------------------------------------

bsg_UCSC_hg38 <- BSgenome.Hsapiens.UCSC.hg38
seq_len_hg38_ucsc=seqlengths(bsg_UCSC_hg38)
#seqinfo(bsg_UCSC)
bsg_NCBI_GRCh38=bsg_UCSC_hg38
seqlevelsStyle(bsg_NCBI_GRCh38) <- "NCBI" #1,2,3,4,..,X,Y
seq_len_GRCh38_NCBI=seqlengths(bsg_NCBI_GRCh38)


bsg_UCSC_hg19 <- BSgenome.Hsapiens.UCSC.hg19
seq_len_hg19_ucsc=seqlengths(bsg_UCSC_hg19)
bsg_NCBI_GRCh37=bsg_UCSC_hg19
seqlevelsStyle(bsg_NCBI_GRCh37) <- "NCBI" #1,2,3,4,..,X,Y
seq_len_GRCh37_NCBI=seqlengths(bsg_NCBI_GRCh37)



## Query GRCh38 from Ensemble --------------------------------------------------------------------
# 

httr_config <- switch(Sys.info()["sysname"],
                      "Linux" = config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"),
                      config())
set_config(config = httr_config)

hub <- AnnotationHub()



## Query Human gene annotations --------------------------------------------

query(hub, c("Homo sapiens", "GRCh38", "Ensembl"))
#edb <- hub[["AH113665"]] #Ensembl 110/111 EnsDb for Homo sapiens ??   

gtf_GRCh38 <- hub[["AH110869"]] #Homo_sapiens.GRCh38.109.gtf 

seqlevels(gtf_GRCh38, pruning.mode="coarse") <- seqlevels(bsg_NCBI_GRCh38)
seqinfo(gtf_GRCh38)=seqinfo(bsg_NCBI_GRCh38)
genome(gtf_GRCh38) <-  "GRCh38"

saveRDS(gtf_GRCh38, file=paste0(data_path,"/gtf_GRCh38_ensembl_109.Rds"))

#gtf_GRCh38=readRDS(file=paste0(data_path,"/gtf_GRCh38_ensembl_109.Rds")) #This is the gtf used

