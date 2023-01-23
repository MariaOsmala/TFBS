library(rtracklayer)

#Database stuff
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
library("dbplyr")
library("dplyr")
#dbplyr is designed to work with database tables as if they were local data frames.
#To demonstrate this I???ll first create an in-memory SQLite database and copy over a dataset:

# Create connection to a database management system

# SQLite only needs a path to the database. (Here, ":memory:" is a special
# path that creates an in-memory database.) Other database drivers
# will require more details (like user, password, host, port, etc.)

GR_list_representative <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches.Rds")


#con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
con <- DBI::dbConnect(RSQLite::SQLite(), "../SQLite/motif-hits.sqlite")

#or should one do this in a loop

df=as.data.frame(GR_list_representative[[1]])
df$PWM=names(GR_list_representative[1])


DBI::dbWriteTable(con, "representative", df )

for(i in 2:length(GR_list_representative)){
  df=as.data.frame(GR_list_representative[[i]])
  df$PWM=names(GR_list_representative[i])
  DBI::dbAppendTable(con, "representative", df)
}

#library(stringr)

representative <- tbl(con, "representative")
#extract data for some motif
# lazily generates query
motif <- representative %>%
  filter(PWM %in% "HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_AGGYGYGANNNNNNNNNNNNNNNTCRTWAA_2_3_YES")

data=DBI::dbSendQuery(con, "SELECT * FROM `representative` WHERE (`PWM` IN ('HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_AGGYGYGANNNNNNNNNNNNNNNTCRTWAA_2_3_YES'))")

#select only seqnames, start, end from the table
res=DBI::dbSendQuery(con, paste0("SELECT seqnames, start, end FROM `representative` WHERE (`PWM` IN ('",names(GR_list_representative[500]),"'))"))
data=DBI::dbFetch(res)
           
DBI::dbClearResult(res)


           names(GR_list_representative[500]) ) 
    
%in% c("UA")

    PWM==names(GR_list_representative[500])) 
str_detect(rowname, "Merc")
# see query
motif %>% show_query()
#> <SQL>
#> SELECT `cyl`, AVG(`mpg`) AS `mpg`
#> FROM `mtcars`
#> GROUP BY `cyl`
#> ORDER BY `mpg` DESC

# execute query and retrieve results
data= motif %>% collect()
#> # A tibble: 3 ?? 2
#>     cyl   mpg
#>   <dbl> <dbl>
#> 1     4  26.7
#> 2     6  19.7
#> 3     8  15.1


#The SQLite database is of size 34.5 GB
# GRangesList object is of size 3.3 in Rds
# 22.3 GB when loaded into memory


DBI::dbDisconnect(con)

DBI::dbWriteTable(con, "mtcars", mtcars)

copy_to(con, mtcars)


DBI::dbWriteTable(con, "mtcars", mtcars)
DBI::dbWriteTable(con, "iris", iris)
DBI::dbListTables(con)

DBI::dbRemoveTable(con, "mtcars")
DBI::dbWriteTable(con, "mtcars", mtcars[1:floor((nrow(mtcars)/2)),])
DBI::dbAppendTable(con, "mtcars", mtcars[(floor((nrow(mtcars)/2))+1):nrow(mtcars),])


#Note that you don???t actually need to load dbplyr with library(dbplyr);
#dplyr automatically loads it for you when it sees you working with a database.
#Database connections are coordinated by the DBI package. Learn more at https://dbi.r-dbi.org/

#Now you can retrieve a table using tbl() (see ?tbl_dbi for more details).
#Printing it just retrieves the first few rows:

mtcars2 <- tbl(con, "mtcars")

#All dplyr calls are evaluated lazily,
#generating SQL that is only sent to the database when you request the data.

# lazily generates query
summary <- mtcars2 %>%
  group_by(cyl) %>%
  summarise(mpg = mean(mpg, na.rm = TRUE)) %>%
  arrange(desc(mpg))

# see query
summary %>% show_query()
#> <SQL>
#> SELECT `cyl`, AVG(`mpg`) AS `mpg`
#> FROM `mtcars`
#> GROUP BY `cyl`
#> ORDER BY `mpg` DESC

# execute query and retrieve results
summary %>% collect()
#> # A tibble: 3 ?? 2
#>     cyl   mpg
#>   <dbl> <dbl>
#> 1     4  26.7
#> 2     6  19.7
#> 3     8  15.1

DBI::dbDisconnect(con)

# Data from Zhang2021 A single cell atlas of chromatin accessibility in the human genome
# * 30 adult and 15 fetal tissues, scATAC-seq
# * 1,152,611 cCREs across 111(adult)+111(fetal)=222 cell types
#   * 890,130 cCREs in adult cell types
#    * 802,025 cCREs in fetal cell types

# cCRE_hg38.tsv.gz
# ================
#
#   The list of 1,154,611 cCREs. There are 7 columns in this file:
#
# 1. chromosome: name of the chromosome.
# 2. hg38_Start: 0-based starting location of the cCRE in the genome (hg38).
# 3. hg38_End: End of the cCRE in the genome (hg38).
# 4. class: Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
# 5. Present in fetal tissues: if this cCRE is detected in at least one fetal tissue
# 6. Present in adult tissues: if this cCRE is detected in at least one adult tissue
# 7. CRE module: The ID of CRE module that the cCRE belongs to.

#What is the definition of cCRE and CRE module
cCRE_hg38=read.table("../CATLAS/cCRE_hg38.tsv.gz",sep="\t", header=FALSE)

Cell_ontology=read.table("../CATLAS/Cell_ontology.tsv",sep="\t", header=TRUE)

#Cell_metadata.tsv.gz
# ====================
#   
#   This file contains metadata for each cell/barcode. It has 6 columns:
#   
# 1. cellID: unique ID.
# 2. logUMI: log10 number of unique fragment passing QC.
# 3. tsse: TSS enrichment
# 4. tissue: The ID of tissue sample
# 5. cell type: The name of the cell type
# 6. Life stage: Fetal or Adult

Cell_metadata=read.table("../CATLAS/Cell_metadata.tsv.gz",sep="\t", header=TRUE)

Celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t", header=TRUE)

# * 435,142 adult cCREs with restricted accessibility in one or a few cell types

#From Mendeley database
files=dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/")

ct_names=do.call(rbind,sapply(files, strsplit, ".bed.gz"))
dimnames(ct_names)=NULL

cCREs_list<-GRangesList()

for(ct in ct_names){
  print(ct)
  cCREs_list[[ct]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/",ct,".bed.gz"),format="bed")


}

saveRDS(cCREs_list,  file = paste0( "../RData/cCREs_list.Rds")) #


sum(sapply(cCREs_list, length)) #1,091,805

#unique
length(unique(unlist(cCREs_list))) #435142

#Hypergeometric test/Fisher's exact test

#A list of 300 000 regions from motif matching. S, m=|S|


#A set of cCREs restrictive to a given cell type. G_0, n=|G_0|

#G: all cCREs across all cell types G, N=|G|

#k: number of cCREs that have motif match

#Which are representative motifs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)

test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])


rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and
#_pdf_spacing_new

#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))
#1062
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject")
GR_list<-readRDS(file = paste0( "../RData/all_motif_matches_sbatch.Rds")) #
length(GR_list)
#3982

GR_list_representative=GR_list[ rep_motifs]
GR_list_representative=GR_list[ which ( ( names(GR_list) %in% rep_motifs)==TRUE)]



#saveRDS(GR_list_representative,  file = paste0( "../RData/representative_motif_matches.Rds")) #


GR_list_representative <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches.Rds")

cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

p_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                dimnames = list(names(GR_list_representative),
                                names(cCREs_list) )) #1062 x 111

e_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
                dimnames = list(names(GR_list_representative),
                                names(cCREs_list) )) #1062 x 111

N=length(unique(unlist(cCREs_list))) #all possible cCREs


## special overlap types
query <- IRanges(c(1, 5, 3, 4), width=c(2, 2, 4, 6))
subject <- IRanges(c(1, 3, 5, 6), width=c(4, 4, 5, 4))

findOverlaps(query, subject, type="start")
findOverlaps(query, subject, type="start", maxgap=1L)
findOverlaps(query, subject, type="end", select="first")
ov <- findOverlaps(query, subject, type="within", maxgap=1L)
ov


for(hits in names(GR_list_representative)){
    print(hits)
    #hits=names(GR_list_representative)[1]
    m=length(GR_list_representative[[hits]])
    for(cell_type in names(cCREs_list)){
        print(cell_type)
        #cell_type=names(cCREs_list)[1]
        #findOverlaps(query=hits, subjects=cCREs_list[cell_type])
      
        #If type is within, the query interval must be wholly contained within 
        #the subject interval. Note that all matches must additionally satisfy 
        #the minoverlap (default 0) constraint described above.
      
        #The maxgap parameter has special meaning with the special overlap types. 
        #For within, it is the maximum amount by which the subject may be wider than the query. 
        #If maxgap is set to -1 (the default), it's replaced internally by 0.
      
      
        fol=findOverlaps(GR_list_representative[[hits]], cCREs_list[[cell_type]], ignore.strand=TRUE)
        #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
        # number of cCREs that have motif match
    
        k=length(unique(fol@to))
    
        #m	the number of white balls in the urn.
        # number of motif matches, 300 000
        
    
    
    
        # N-m
        # the number of black balls in the urn.
        # the number of all possible cCREs - the number of motif matches
    
        Nm=N-m
    
        # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
        # The number of cell type specific cCREs
    
        n=length(cCREs_list[[cell_type]])
    
        e = log2( (k / n) / (m / N) )
    
        #k=oFg_hs - 1, m=oBg_hs, N-m=nBg_hs - oBg_hs n= nFg_hs
      
        #oFg = fromIntegral $ length $ filter (`S.member` motifs) fg :: Double
        #oBg = fromIntegral $ S.size motifs :: Double
        #nFg = fromIntegral $ length fg :: Double
        #nBg = fromIntegral nBg' :: Double
      
        if(e >= 0){
          p <- phyper(k - 1, m, Nm, n, lower.tail=F, log.p=T) / log(10)
    
        }else{
          # lower.tail=TRUE (default), probabilities are P[X<=k], otherwise, P[X > k].
          #k in 0, ..., n
          # Why the lower tail TRUE? Should be false?
          p <- phyper(k, m, Nm, n, lower.tail=T, log.p=T) / log(10)
        }
        
        
        p_matrix[hits, cell_type]=p
        e_matrix[hits, cell_type]=3
    
  
  }
  
}
  
saveRDS(p_matrix,  file = paste0( "../RData/p_matrix.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/e_matrix.Rds")) #
  
  
  
#Fisher's exact test for association:
#Are the motif hits for a given TF found in the cell-type specific cCREs more often than would be expected by change?
#NUll hypothesis: cell-type specific cCREs are independent of the motif matches, the number of overlaps follows a hypergeometric distribution

p(overlap=k)=

?phyper

phyper(k, m, Nm, n, lower.tail = TRUE, log.p = FALSE)
k
m
Nm
n
# k	vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
# number of cCREs that have motif match

#m	the number of white balls in the urn.
# number of motif matches, 300 000

# N-m
# the number of black balls in the urn.
# the number of all possible cCREs - the number of motif matches

# n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
# The number of cell type specific cCREs

#p	probability, it must be between 0 and 1.



#How the enrichment analysis performed
# enrichment :: S.HashSet B.ByteString -> [B.ByteString] -> Int -> IO (Double, Double)
# enrichment motifs fg nBg' = if e >= 0
#     then R.runRegion $ do
#         p <- [r| phyper(oFg_hs - 1, oBg_hs, nBg_hs - oBg_hs, nFg_hs, lower.tail=F, log.p=T) / log(10)
#             |]
#         R.unSomeSEXP p $ \x -> case H.hexp x of
#             H.Real x'  -> return (e, head $ toList x')
#     else R.runRegion $ do
#         p <- [r| phyper(oFg_hs, oBg_hs, nBg_hs - oBg_hs, nFg_hs, lower.tail=T, log.p=T) / log(10)
#             |]
#         R.unSomeSEXP p $ \x -> case H.hexp x of
#             H.Real x'  -> return (e, head $ toList x')
#   where
#     e = logBase 2 $ (oFg / nFg) / (oBg / nBg)
#     oFg = fromIntegral $ length $ filter (`S.member` motifs) fg :: Double
#     oBg = fromIntegral $ S.size motifs :: Double
#     nFg = fromIntegral $ length fg :: Double
#     nBg = fromIntegral nBg' :: Double



# * Phylogenetic analysis to group fetal and adult cell types into different cell lineages, enrichment for lineage-specific TFBSs?
#   * Differentially accessible (DA) elements between adults and fetal cell lines/cell lineage groups.
#    Adult and fetal-specific cCREs for the cell-lineage groups. Specific cCREs underlying fetal or adult-specific regulatory programs.

#WHERE IS THIS DATA?

cCRE_hg38=read.table("../CATLAS/cCRE_hg38.tsv.gz",sep="\t", header=FALSE)
nrow(cCRE_hg38)

Cell_ontology=read.table("../CATLAS/Cell_ontology.tsv",sep="\t", header=TRUE)



# 1. chromosome: name of the chromosome.
# 2. hg38_Start: 0-based starting location of the cCRE in the genome (hg38).
# 3. hg38_End: End of the cCRE in the genome (hg38).
# 4. class: Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
# 5. Present in fetal tissues: if this cCRE is detected in at least one fetal tissue
# 6. Present in adult tissues: if this cCRE is detected in at least one adult tissue
# 7. CRE module: The ID of CRE module that the cCRE belongs to.

#1,152,611 cCREs across 111(adult)+111(fetal)=222 cell types
library(readr)
gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")

cCREs=read_delim("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_cCRE_list_GRCh38.tsv.gz", delim="\t", col_names=TRUE)
nrow(cCREs) #1154611, this is the same table as cCRE_hg38

names(cCREs)=c("chr", "start", "end", "class", "fetal", "adult", "CREmoduleID")

library("tidyverse")
library("GenomicRanges")
modules=cCREs %>%
  group("CREmoduleID")

by_CREdmodules <- cCREs %>% group_by(CREmoduleID)

by_CREdmodules %>% tally(sort=TRUE)

group_indices=by_CREdmodules %>% group_indices()

CREmodules <- GRangesList()

for(gk in 1:150){
  print(gk)
  
  data=as.data.frame(cCREs %>% filter( CREmoduleID==gk) %>% select(c("chr", "start", "end")))
  tmp_GRanges=GRanges(seqnames = Rle( data$chr, rep(1, nrow(data)) ),
                      ranges = IRanges(start=data$start, end = data$end) ,
                      strand = Rle(strand( rep("*", nrow(data))) , rep(1, nrow(data))  ))
  
  
  isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  
  seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  
  #try dropSeqlevels or ?keepStandardChromosomes from GenomicRanges library.
  
  
  tmp_GRanges=dropSeqlevels(tmp_GRanges, names(seqlengths(tmp_GRanges)[ which( (as.vector(names(seqlengths(tmp_GRanges))) %in% names(seqlengths)==FALSE) ) ]),
                            pruning.mode = "coarse")
  
  #tmp_GRanges=tmp_GRanges[which( (as.vector(seqnames(tmp_GRanges)) %in% names(seqlengths)==FALSE) )]
  
  #seqlengths(tmp_GRanges)= seqlengths(tmp_GRanges)[ which( (as.vector(names(seqlengths(tmp_GRanges))) %in% names(seqlengths)==FALSE) ) ]
  
  #CREmodules may contain some chromosomes not in gtf
  
  
  seqlengths(tmp_GRanges)=seqlengths[names( seqlengths(tmp_GRanges))]
  
  genome(tmp_GRanges) <- "GRCh38"
  
  CREmodules[[gk]]=tmp_GRanges
  
  
}


#saveRDS(CREmodules,  file ="../RData/CREmodules.Rds") 
CREmodules<-readRDS(file = "../RData/CREmodules.Rds")

#Modules connected to cell_types, groups of cell types

Meta_modules=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5D_CRE_Module_Motif_enrichment.csv",sep="\t")
Meta_modules=as.data.frame(t(Meta_modules[c(1,2), -1]))
rownames(Meta_modules)=NULL
names(Meta_modules)=c("ModuleGroup", "ModuleID")
Meta_modules[,"ModuleID"]=as.numeric(Meta_modules[,"ModuleID"])

saveRDS(Meta_modules,  file ="../RData/Meta_modules.Rds") 

qn=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Quantile_normalized_log_RPKM.tsv-NN6XGL.gz",sep="\t", header=TRUE)

RAW_RPKM=read.table("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Raw_RPKM.tsv.gz",sep="\t", header=TRUE)

cell_types=colnames(RAW_RPKM)[-1]

#RAW_RPKM=read_delim("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 5/5A_Raw_RPKM.tsv.gz", delim="\t", col_names=TRUE)
RAW_RPKM[1:5,]

# * 435,142 adult cCREs with restricted accessibility in one or a few cell types

#From Mendeley database
files=dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/")

ct_names=do.call(rbind,sapply(files, strsplit, ".bed.gz"))
dimnames(ct_names)=NULL

cCREs_list<-GRangesList()

for(ct in ct_names){
  print(ct)
  cCREs_list[[ct]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/",ct,".bed.gz"),format="bed")
  
  
}


# * Cell-type specificity of cCREs across the cell types.
# 1.2 million cCREs clustered into 150 cis-regulatory modules (CRMs).
# Putative master regulators of fetal and adult cell types.



#Aims:

#Motif enrichment analysis on cell-type restricted cCREs
  #345,142 cCREs that demonstrated restricted accessibility in one or a few (how many?) cell types
  #to reveal putative TFs of each cell type
  #cCREs restricted to a certain cell type, 10 most significant TFs selected of a single cell type
  #results is cell type x TF matrix
  # Color represents -log_10(P-val), depleted (blue), enriched (red)


#cCRE_by_cell_type
#=================
#This directory stores cCRE by cell type, containing three files:

#1. matrix.mtx.gz: matrix in matrix market format, each line represents a cCRE-cell type pair. For example, "2 22" means cCRE No. 2 is accessibile in cell type No. 22.
#2. celltypes.txt.gz: column names (cell types) of the matrix.
#3. cCREs.bed.gz: row names (cCREs) of the matrix.
library("Matrix")
matrix=readMM("../CATLAS/matrix.tsv.gz") #this is 1,154,611 x 222 matrix, contains TRUE or FALSE

tmp=do.call(rbind, apply(matrix,1,table)) #For each cCREs, in how many it is present
restricted_cCREs=which(tmp[,"TRUE"]==1) #cCREs that are present in a single cell type  347991 (paper: 345142 cCREs, maybe blacklist removed)

celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t")

#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)



cCREs=import("../CATLAS/cCREs.bed.gz",format="bed")





#read TF peaks
#Yimengs motifs
#RFX3 GLI3
#FLI1 GLI3
#TCF7

TCF7=read.table("../../Results/MOODS/TCF7.csv",sep=",")

#computed the odds ratio and the significance of enrichment in each cluster,
#comparing to a non-dynamic bin pool using Fisher???s exact test.
#The non-dynamic bin pool was sampled with replacement to
#match the distribution of average signal strength from the dynamic bins.
#Following that, significant TF PWMs were grouped in subfamilies using
#the structural information from TFClass71 because they share similar
#if not identical binding motifs. The top significantly over-represented TFs
#and their associated subfamilies were reported.


#Motif enrichment and GO analysis. Motif enrichment for each cell type.
# Motif enrichment for each cell type and histone modifications
# were carried out using ChromVAR61. Briefly, mapped reads were converted to
# cell-to-bin matrices with a bin size of 1,000 bp for four histone profiles.
# Reads for each bin were summarized from all cells of the same
# groups from transcriptome-based clustering. GC bias and background peaks were
# calculated and motif enrichment score for each cell type was
# then computed using the computeDeviations function of ChromVAR61.
# Motif enrichment for each CRE module: Motif enrichment for each CRE module
# was analyzed using Homer (v.4.11)62. We scanned a region of ??200 bp
# around the center of the element for both de novo and known motif enrichment analysis
# (from the JASPAR database63). We used the total peak list as the background
# for motif enrichment analysis of cCREs in each group.
#
#61.Schep, A. N., Wu, B., Buenrostro, J. D. & Greenleaf, W. J. chromVAR: inferring
#transcription-factor-associated accessibility from single-cell epigenomic data. Nat. Methods 14, 975???978 (2017).

#62.Heinz, S. et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory
#elements required for macrophage and B cell identities. Mol. Cell 38, 576???589 (2010).
