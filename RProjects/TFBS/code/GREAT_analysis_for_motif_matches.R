library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

.libPaths("/projappl/project_2006203/project_rpackages_4.3.0" )

data_path="/scratch/project_2006203/TFBS/"
#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]



length(unique(rep_motifs))

representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")


motifs=names(representative_motif_matches)


grep("HNF4", motifs)

library("rGREAT") #already installed in csc

for(i in 1:length(motifs)){
  print(i)
  test=representative_motif_matches[[i]]
  strand(test)="*"
  genome(test)="hg38"
  
  res2=great(test, #input regions, chromosome names should be consistent between these and TSS, try getTSS
        gene_sets="GO:BP", #default supported gene sets collections, see Genesets
        tss_source="hg38", #"TxDb.Hsapiens.UCSC.hg38.knownGene", #see TSS section
        #biomart_dataset = NULL, #BioMartGOGeneSets::supportedOrganisms.
        min_gene_set_size = 5, 
        mode = "basalPlusExt", #the mode to extent genes, basalPlusExt', 'twoClosest' and 'oneClosest' See extendTSS
        basal_upstream = 5000, #Only in basalPlusExt mode, number of base pairs extending to the upstream of TSS to form the basal domains.
        basal_downstream = 1000, #Only in basalPlusExt mode, number of base pairs extending to the downstream of TSS to form the basal domains.
        extension = 1000000, #Extensions from the basal domains.
        #extended_tss = NULL, 
        #background = NULL, 
        exclude = "gap", #gap regions (which can be get by getGapFromUCSC).
        cores = 1, #number of cores to use
        verbose = great_opt$verbose)
  
  tb = getEnrichmentTable(res2)
  head(tb)
  
  plotVolcano(res2)
  plotRegionGeneAssociations(res2)
  getRegionGeneAssociations(res2)
  #plotRegionGeneAssociations(res, term_id = "HALLMARK_APOPTOSIS")
  #getRegionGeneAssociations(res, term_id = "HALLMARK_APOPTOSIS")
  shinyReport(res2)
  
  
  
}


job = submitGreatJob(gr)
job = submitGreatJob(gr, species = "mm9") # of course, gr should be from mm9
job = submitGreatJob(gr, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
job = submitGreatJob(gr, rule = "twoClosest", adv_twoDistance = 2000)
job = submitGreatJob(gr, rule = "oneClosest", adv_oneDistance = 2000)

#Also you can choose different versions of GREAT for the analysis.

job = submitGreatJob(gr, version = "3.0")
job = submitGreatJob(gr, version = "2.0")