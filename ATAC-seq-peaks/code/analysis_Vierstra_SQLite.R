library("rtracklayer")
library("dbplyr")
library("dplyr")

#Database stuff
setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

gtf<-readRDS("../../RProjects/TFBS/gtf.Rds")

#How the database was created, Vierstra TFBS/RProjects/TFBS/MOODS_results_to_database.R


con <- DBI::dbConnect(RSQLite::SQLite(), "../../Results/MOODS_Vierstra_SQLite/motif-hits.sqlite")



#motif.hits <- tbl(con, "motif-hits")

motifs <- as_tibble(tbl(con, "motifs")) #3980, should be 3982, 2 missing?

#Which motifs are representative

#Which are representative motifs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))
#1062

#Representative motif ids

#rep_id= motifs %>% filter("motif_id", PWMs %in% rep_motifs)

#data=as.data.frame(motifs %>% filter( PWMs %in% rep_motifs ) %>% select("motif_id" ) )

rep_id=as_tibble(motifs %>% filter( PWMs %in% rep_motifs ) %>% select(c("PWMs","motif_id")) )

#for which motifs the query failed

rep_motifs[which((rep_motifs %in% rep_id$PWMs)==FALSE)]

length(which((rep_motifs %in% rep_id$PWMs)==FALSE)) #0


cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

p_matrix=matrix(nrow=length(rep_id$PWMs), ncol=length(cCREs_list),
                dimnames = list( rep_id$PWMs,
                                names(cCREs_list) )) #1062 x 111

e_matrix=matrix(nrow=length(rep_id$PWMs), ncol=length(cCREs_list),
                dimnames = list(rep_id$PWMs,
                                names(cCREs_list) )) #1062 x 111

N=length(unique(unlist(cCREs_list))) #all possible cCREs #435142

#m: a list of all cCREs that have a motif match



#

#select only seqnames, start, end from the table
#head(data.frame(rep_motifs, sort(as.vector(rep_id$PWMs)) ))
#table(rep_motifs==sort(as.vector(rep_id$PWMs)) ) all TRUE

representatives<-GRangesList()


for(hits in as.vector(rep_id$PWMs)){
  print(hits)
  #hits=as.vector(rep_id$PWMs)[1]

  #id of the motif
  id=rep_id %>% filter(PWMs==hits) %>% select("motif_id")
  
  res=DBI::dbSendQuery(con, paste0("SELECT chr, position, sequence, strand, score FROM `motif-hits` WHERE (`motif_id` IN ('", id ,"'))"))
  
  data=DBI::dbFetch(res)
  
  DBI::dbClearResult(res)
  
  #   tmp=tb[tb$PWM==pwm,]
  tmp_GRanges=GRanges(seqnames = Rle( data$chr, rep(1, nrow(data)) ),
     ranges = IRanges(start=data$position+1, end = data$position+ nchar(data$sequence[1]) ),
     strand = Rle(strand(data$strand  ), rep(1, nrow(data))  ), score=data$score)
  
  

  # 
  isCircular(seqinfo(tmp_GRanges))=rep(FALSE, length(seqinfo(tmp_GRanges)))
  # 
  seqlengths=seqlengths(seqinfo(gtf))[c(as.character(seq(1,22,1)), "X", "Y")]
  names(seqlengths)=paste0("chr", names(seqlengths))
  # 
  seqlengths(tmp_GRanges)=seqlengths
  genome(tmp_GRanges) <- "GRCh38"
  # 
  #   #choose 300,000 with the best score
  # 
  tmp_GRanges=tmp_GRanges[order(tmp_GRanges$score, decreasing=TRUE),]
  # 
  representatives[[ hits ]]=tmp_GRanges
  
  m=length(tmp_GRanges) #138061
  
  #unlist(lapply(cCREs_list, length))
  
  #library(ggplot2)
  # Basic histogram
  
  #df=data.frame(cCREs_nro= unlist(lapply(cCREs_list, length)))
  
  #max(df$cCREs_nro) #48391
  #min(df$cCREs_nro) #846
  
  #ggplot(df, aes(x=_nro)) + geom_histogram()
  # Change the width of bins
  #pdf(file = "../Figures/cCREs-nros.pdf" )
  #width, height, onefile, family, title, fonts, version,
  #paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
  #useDingbats, useKerning, fillOddEven, compress)
  #ggplot(df, aes(x=cCREs_nro)) + geom_histogram(binwidth=1000)
  #dev.off()
  
  #Width of cCREs, these are all the same 400 bp?
  
  #widths=lapply(cCREs_list, width)
  #lapply(widths, table)
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #unlist(lapply(cCREs_list, length))
    
    #findOverlaps(query, subject)
    fol=findOverlaps(tmp_GRanges, cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
    #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
    # number of cCREs that have motif match
    
    k=length(unique(fol@to)) #This is now number of cCREs that have at least one motif hit #245
    #k=length(unique(fol@from)) This is now number of hits that overlap with cCREs (there can be multiple motifs hits in a single CRE)
    
    #m	the number of white balls in the urn.
    # number of motif matches, 300 000
    
    # N-m
    # the number of black balls in the urn.
    # the number of all possible cCREs - the number of motif matches
    
    Nm=N-m #list of all CREs that do not have a motif
    
    # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
    # The number of cell type specific cCREs
    
    n=length(cCREs_list[[cell_type]]) #This is correct
    
    #Hypergeometric test/Fisher's exact test
    
    #m: a list of all cCREs that have a motif match
    
    #n: A set of cCREs restrictive to a given cell type. G_0, n=|G_0|
    
    #N: G: all cCREs across all cell types G, N=|G|
    
    #k: number of cell-line-specific cCREs that have motif match
    
    
    fe=FishersExact(k,n,m,N)
    p_matrix[hits, cell_type]=fe$p
    e_matrix[hits, cell_type]=fe$e
    
    
    
  }
  
}


DBI::dbDisconnect(con)

saveRDS(representatives,  file = paste0( "../RData/Vierstra_representatives.Rds")) #
saveRDS(p_matrix,  file = paste0( "../RData/p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/e_matrix_human_cell_type_restricted.Rds")) #







#motif %>% show_query()
#> <SQL>
#> SELECT `cyl`, AVG(`mpg`) AS `mpg`
#> FROM `mtcars`
#> GROUP BY `cyl`
#> ORDER BY `mpg` DESC

# execute query and retrieve results
#data= motif %>% collect()


#or should one do this in a loop



#library(stringr)

#representative <- tbl(con, "representative")
#extract data for some motif
# lazily generates query
#motif <- representative %>%
#  filter(PWM %in% "HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_AGGYGYGANNNNNNNNNNNNNNNTCRTWAA_2_3_YES")

#data=DBI::dbSendQuery(con, "SELECT * FROM `representative` WHERE (`PWM` IN ('HOXD12_EOMES_CAP-SELEX_TGAGTC40NAAC_AAB_AGGYGYGANNNNNNNNNNNNNNNTCRTWAA_2_3_YES'))")

#select only seqnames, start, end from the table
# res=DBI::dbSendQuery(con, paste0("SELECT seqnames, start, end FROM `representative` WHERE (`PWM` IN ('",names(GR_list_representative[500]),"'))"))
# data=DBI::dbFetch(res)
#            
# DBI::dbClearResult(res)


           # names(GR_list_representative[500]) ) 
    
# %in% c("UA")
# 
#     PWM==names(GR_list_representative[500])) 
# str_detect(rowname, "Merc")
# see query
# motif %>% show_query()
#> <SQL>
#> SELECT `cyl`, AVG(`mpg`) AS `mpg`
#> FROM `mtcars`
#> GROUP BY `cyl`
#> ORDER BY `mpg` DESC

# execute query and retrieve results
# data= motif %>% collect()
#> # A tibble: 3 ?? 2
#>     cyl   mpg
#>   <dbl> <dbl>
#> 1     4  26.7
#> 2     6  19.7
#> 3     8  15.1


#The SQLite database is of size 34.5 GB
# GRangesList object is of size 3.3 in Rds
# 22.3 GB when loaded into memory


# DBI::dbDisconnect(con)
# 
# DBI::dbWriteTable(con, "mtcars", mtcars)
# 
# copy_to(con, mtcars)
# 
# 
# DBI::dbWriteTable(con, "mtcars", mtcars)
# DBI::dbWriteTable(con, "iris", iris)
# DBI::dbListTables(con)
# 
# DBI::dbRemoveTable(con, "mtcars")
# DBI::dbWriteTable(con, "mtcars", mtcars[1:floor((nrow(mtcars)/2)),])
# DBI::dbAppendTable(con, "mtcars", mtcars[(floor((nrow(mtcars)/2))+1):nrow(mtcars),])


#Note that you don???t actually need to load dbplyr with library(dbplyr);
#dplyr automatically loads it for you when it sees you working with a database.
#Database connections are coordinated by the DBI package. Learn more at https://dbi.r-dbi.org/

#Now you can retrieve a table using tbl() (see ?tbl_dbi for more details).
#Printing it just retrieves the first few rows:

#mtcars2 <- tbl(con, "mtcars")

#All dplyr calls are evaluated lazily,
#generating SQL that is only sent to the database when you request the data.

# lazily generates query
# summary <- mtcars2 %>%
#   group_by(cyl) %>%
#   summarise(mpg = mean(mpg, na.rm = TRUE)) %>%
#   arrange(desc(mpg))
# 
# see query
# summary %>% show_query()
#> <SQL>
#> SELECT `cyl`, AVG(`mpg`) AS `mpg`
#> FROM `mtcars`
#> GROUP BY `cyl`
#> ORDER BY `mpg` DESC

# execute query and retrieve results
# summary %>% collect()
#> # A tibble: 3 ?? 2
#>     cyl   mpg
#>   <dbl> <dbl>
#> 1     4  26.7
#> 2     6  19.7
#> 3     8  15.1

# DBI::dbDisconnect(con)

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
# cCRE_hg38=read.table("../CATLAS/cCRE_hg38.tsv.gz",sep="\t", header=FALSE)
# 
# Cell_ontology=read.table("../CATLAS/Cell_ontology.tsv",sep="\t", header=TRUE)



# * 435,142 adult cCREs with restricted accessibility in one or a few cell types

#From Mendeley database
# files=dir("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/")

# ct_names=do.call(rbind,sapply(files, strsplit, ".bed.gz"))
# dimnames(ct_names)=NULL

# cCREs_list<-GRangesList()

# for(ct in ct_names){
#   print(ct)
#   cCREs_list[[ct]]=import(paste0("../CATLAS/yv4fzv6cnm-4/Zhang et al Figure 2/2E_Cell_type_restricted_peaks/",ct,".bed.gz"),format="bed")
# 
# 
# }
# 
# saveRDS(cCREs_list,  file = paste0( "../RData/cCREs_list.Rds")) #
# 
# 
# sum(sapply(cCREs_list, length)) #1,091,805
# 
#unique
# length(unique(unlist(cCREs_list))) #435142

#Hypergeometric test/Fisher's exact test

#A list of 300 000 regions from motif matching. S, m=|S|


#A set of cCREs restrictive to a given cell type. G_0, n=|G_0|

#G: all cCREs across all cell types G, N=|G|

#k: number of cCREs that have motif match

#Which are representative motifs

# representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
# 
# test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
# 
# 
# rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and
#_pdf_spacing_new

#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

# length(unique(rep_motifs))
#1062
# setwd("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RProject")
# GR_list<-readRDS(file = paste0( "../RData/all_motif_matches_sbatch.Rds")) #
# length(GR_list)
#3982

# GR_list_representative=GR_list[ rep_motifs]
# GR_list_representative=GR_list[ which ( ( names(GR_list) %in% rep_motifs)==TRUE)]



#saveRDS(GR_list_representative,  file = paste0( "../RData/representative_motif_matches.Rds")) #


# GR_list_representative <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches.Rds")

# cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
# length(cCREs_list)
# 
# p_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
#                 dimnames = list(names(GR_list_representative),
#                                 names(cCREs_list) )) #1062 x 111
# 
# e_matrix=matrix(nrow=length(GR_list_representative), ncol=length(cCREs_list),
#                 dimnames = list(names(GR_list_representative),
#                                 names(cCREs_list) )) #1062 x 111
# 
# N=length(unique(unlist(cCREs_list))) #all possible cCREs
# 

## special overlap types
# query <- IRanges(c(1, 5, 3, 4), width=c(2, 2, 4, 6))
# subject <- IRanges(c(1, 3, 5, 6), width=c(4, 4, 5, 4))
# 
# findOverlaps(query, subject, type="start")
# findOverlaps(query, subject, type="start", maxgap=1L)
# findOverlaps(query, subject, type="end", select="first")
# ov <- findOverlaps(query, subject, type="within", maxgap=1L)
# ov
# 

# for(hits in names(GR_list_representative)){
#     print(hits)
#     #hits=names(GR_list_representative)[1]
#     m=length(GR_list_representative[[hits]])
#     for(cell_type in names(cCREs_list)){
#         print(cell_type)
        #cell_type=names(cCREs_list)[1]
        #findOverlaps(query=hits, subjects=cCREs_list[cell_type])
      
        #If type is within, the query interval must be wholly contained within 
        #the subject interval. Note that all matches must additionally satisfy 
        #the minoverlap (default 0) constraint described above.
      
        #The maxgap parameter has special meaning with the special overlap types. 
        #For within, it is the maximum amount by which the subject may be wider than the query. 
        #If maxgap is set to -1 (the default), it's replaced internally by 0.
      
      
        # fol=findOverlaps(GR_list_representative[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
        #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
        # number of cCREs that have motif match
    
        # k=length(unique(fol@to))
    
        #m	the number of white balls in the urn.
        # number of motif matches, 300 000
        
    
    
    
        # N-m
        # the number of black balls in the urn.
        # the number of all possible cCREs - the number of motif matches
    
        # Nm=N-m
    
        # n	the number of balls drawn from the urn, hence must be in 0,1,\dots, m+n0,1,???,m+n=N
        # The number of cell type specific cCREs
    
        # n=length(cCREs_list[[cell_type]])
    
        # e = log2( (k / n) / (m / N) )
    
        #k=oFg_hs - 1, m=oBg_hs, N-m=nBg_hs - oBg_hs n= nFg_hs
      
        #oFg = fromIntegral $ length $ filter (`S.member` motifs) fg :: Double
        #oBg = fromIntegral $ S.size motifs :: Double
        #nFg = fromIntegral $ length fg :: Double
        #nBg = fromIntegral nBg' :: Double
      
        # if(e >= 0){
        #   p <- phyper(k - 1, m, Nm, n, lower.tail=F, log.p=T) / log(10)
        # 
        # }else{
        #   # lower.tail=TRUE (default), probabilities are P[X<=k], otherwise, P[X > k].
        #   #k in 0, ..., n
        #   # Why the lower tail TRUE? Should be false?
        #   p <- phyper(k, m, Nm, n, lower.tail=T, log.p=T) / log(10)
        # }
        # 
        # 
        # p_matrix[hits, cell_type]=p
        # e_matrix[hits, cell_type]=3
        # 
  
  # }
  
# }
  
# saveRDS(p_matrix,  file = paste0( "../RData/p_matrix.Rds")) #
# saveRDS(e_matrix,  file = paste0( "../RData/e_matrix.Rds")) #
  
  
  
#Fisher's exact test for association:
#Are the motif hits for a given TF found in the cell-type specific cCREs more often than would be expected by change?
#NUll hypothesis: cell-type specific cCREs are independent of the motif matches, the number of overlaps follows a hypergeometric distribution

# p(overlap=k)=
# 
# ?phyper
# 
# phyper(k, m, Nm, n, lower.tail = TRUE, log.p = FALSE)
# k
# m
# Nm
# n
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
# 





# * Phylogenetic analysis to group fetal and adult cell types into different cell lineages, enrichment for lineage-specific TFBSs?
#   * Differentially accessible (DA) elements between adults and fetal cell lines/cell lineage groups.
#    Adult and fetal-specific cCREs for the cell-lineage groups. Specific cCREs underlying fetal or adult-specific regulatory programs.




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
# library("Matrix")
# matrix=readMM("../CATLAS/matrix.tsv.gz") #this is 1,154,611 x 222 matrix, contains TRUE or FALSE
# 
# tmp=do.call(rbind, apply(matrix,1,table)) #For each cCREs, in how many it is present
# restricted_cCREs=which(tmp[,"TRUE"]==1) #cCREs that are present in a single cell type  347991 (paper: 345142 cCREs, maybe blacklist removed)
# 
# celltypes=read.table("../CATLAS/celltypes.txt.gz",sep="\t")
# 
# #source("https://bioconductor.org/biocLite.R")
# #biocLite("rtracklayer")
# library(rtracklayer)
# 
# 
# 
# cCREs=import("../CATLAS/cCREs.bed.gz",format="bed")
# 
# 



#read TF peaks
#Yimengs motifs
#RFX3 GLI3
#FLI1 GLI3
#TCF7

# TCF7=read.table("../../Results/MOODS/TCF7.csv",sep=",")

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
