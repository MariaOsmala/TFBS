#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))
#.libPaths("/projappl/project_2006203/project_rpackages_4.2.1", include.site=FALSE)
#.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new", include.site=FALSE)


#install.packages("BiocManager")
#BiocManager::install("rtracklayer") Do not install this to personal library, things break
library("rtracklayer")
#BiocManager::install("TFBSTools")
library("TFBSTools")
#BiocManager::install("Biostrings")
library("Biostrings")
#BiocManager::install("motifmatchr")
library("motifmatchr")
#BiocManager::install("GenomicRanges")
library("GenomicRanges")
library("readr")
library("ggplot2")
#install.packages("ggbreak")
.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new")
library(ggbreak)

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# wget https://packagemanager.posit.co/bioconductor/packages/3.17/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz
# wget https://packagemanager.posit.co/bioconductor/packages/3.17/data/annotation/src/contrib/BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000.tar.gz
# wget https://packagemanager.posit.co/bioconductor/packages/3.17/data/annotation/src/contrib/BSgenome.Mmusculus.UCSC.mm39_1.4.3.tar.gz
# wget https://packagemanager.posit.co/bioconductor/packages/3.17/data/annotation/src/contrib/BSgenome.Mmusculus.UCSC.mm10_1.4.3.tar.gz


# install.packages("/projappl/project_2006472/project_rpackages_4.3.0_new/BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000.tar.gz", repos = NULL, type = "source")
# install.packages("/projappl/project_2006472/project_rpackages_4.3.0_new/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", repos = NULL, type = "source")
# install.packages("/projappl/project_2006472/project_rpackages_4.3.0_new/BSgenome.Mmusculus.UCSC.mm10_1.4.3.tar.gz", repos = NULL, type = "source")
# install.packages("/projappl/project_2006472/project_rpackages_4.3.0_new/BSgenome.Mmusculus.UCSC.mm39_1.4.3.tar.gz", repos = NULL, type = "source")


#library("BSgenome.Hsapiens.UCSC.hg38")
#library("BSgenome.Hsapiens.NCBI.GRCh38")
#library("BSgenome.Mmusculus.UCSC.mm39")
#library("BSgenome.Mmusculus.UCSC.mm10")

library(BSgenome.Mmusculus.UCSC.mm39)
bsg <- BSgenome.Mmusculus.UCSC.mm39
seq_len_mm39=seqlengths(bsg)
names(seq_len_mm39)=gsub("chr", "", names(seq_len_mm39))
seqinfo(bsg)
seqlevelsStyle(bsg) <- "UCSC" #1,2,3,4,..,X,Y
saveRDS(seq_len_mm39, file="/scratch/project_2006472/RData/mm39_chrom_lens.Rds")

library(BSgenome.Mmusculus.UCSC.mm10)
bsg_mm10 <- BSgenome.Mmusculus.UCSC.mm10
seq_len_mm10=seqlengths(bsg_mm10)
names(seq_len_mm10)=gsub("chr", "", names(seq_len_mm10))
seqinfo(bsg_mm10)
seqlevelsStyle(bsg_mm10) <- "UCSC" #1,2,3,4,..,X,Y
saveRDS(seq_len_mm10, file="/scratch/project_2006472/RData/mm10_chrom_lens.Rds")


#hub <- AnnotationHub()
#query(hub, c("Homo sapiens", "GRCh37", "Ensembl"))
#test=query(hub, c("Homo sapiens", "GRCh38", "Ensembl"))
#gtf <- hub[["AH10684"]] #Homo_sapiens.GRCh37.75.gtf
#saveRDS(gtf, file="/scratch/project_2006472/RData/gtf_GRCh37.Rds")


#hub <- AnnotationHub()
#query(hub, c("Homo sapiens", "GRCh37", "Ensembl"))
#test=query(hub, c("Homo sapiens", "GRCh38", "Ensembl"))
#gtfGRCh38 <- hub[["AH110869"]] #Homo_sapiens.GRCh38.109.gtf 
#saveRDS(gtfGRCh38, file="/scratch/project_2006472/RData/gtf_GRCh38.Rds")



#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

#wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToHg38.over.chain.gz
#wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz



data_path="/scratch/project_2006203/TFBS/"
#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives[which(representatives$new_representative=="YES"),]
#length(unique(rep_motifs$ID))#1031

######################### REPRESENTATIVE MOUSE mm39 MOTIF MATCHES ####################################################
# top_motif_matches_mouse_final=readRDS(file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_mouse_final.Rds") # 3124? should be 3294, more?
# match_numbers=unlist(lapply(top_motif_matches_mouse_final, length))
# saveRDS(match_numbers,"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_top_motif_matches_mouse_final.Rds") 
# df=data.frame(ID=names(match_numbers),match_numbers=as.vector(match_numbers))
# merged_representatives <- merge(representatives, df[,c('ID','match_numbers')], by = 'ID', all.x = TRUE)
# floor(min(top_motif_matches_mouse_final[[1]]$score))
# f1 <- function(x) floor(min(x$score))
# thresholds=lapply(top_motif_matches_mouse_final, f1) #, simplity="array"
# df=data.frame(ID=names(top_motif_matches_mouse_final),threshold=unlist(thresholds))
# merged_representatives2 <- merge(merged_representatives, df[,c('ID','threshold')], by = 'ID', all.x = TRUE)
# write.table(merged_representatives2, "../../PWMs_final/metadata_representatives_mouse_match_numbers_thresholds.tsv",sep="\t", 
#              col.names=TRUE, row.names=FALSE) #3294
# pdf(file = paste0(data_path, "ATAC-seq-peaks/Figures/motif-match-number-distribution-mouse-final.pdf" ))
# df=data.frame(match_nro=match_numbers)
# ggplot(df, aes(x=match_nro)) + geom_histogram(color="black", fill="white",aes(y = after_stat(count / sum(count))),binwidth=20000)+
#   scale_y_cut(breaks=c(0.05), which=c(1,2), scales=c(0.5,0.5)) + ylim(c(0, 1))+#xlim(c(280000,325000))
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+   xlim(0, 6e5)+ scale_x_continuous(labels = scales::comma)+
#   labs(x='Number of matches', y='Frequency', title='Histogram of motif match numbers, binwidth 20000')+ theme_minimal()
# dev.off()
# representative_motif_matches=top_motif_matches_mouse_final[rep_motifs]
# saveRDS(representative_motif_matches, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_mouse_final.Rds")

#https://www.ncbi.nlm.nih.gov/gene/105804841 GRCh38 chr7:156,790,708-156,793,079 #this is 2,372 bp.

#coordinates in GRCh37 chr7:156583402..156585773

#Vista browser, this is in hg19:  chr7:156,583,782-156,584,569

gtf_GRCh38 <- readRDS("/scratch/project_2006203/TFBS/RProjects/TFBS/gtf.Rds")

# Get motif matches for example motifs in peaks 

#The VISTA enhancer, 787 bp region
#156584569-156583782 This is hg19, convert to GRCh38
test_GRanges=GRanges(seqnames = Rle( "chr7", rep(1, 1) ),
                     ranges = IRanges(start=156583782, end = 156584569) )

genome(test_GRanges) <- "GRCh37"

ChainFilePath='/projappl/project_2006203/liftOver/hg19ToHg38.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
test_GRanges_GRCh38=unlist(rtracklayer::liftOver(test_GRanges, ch))
isCircular(seqinfo(test_GRanges_GRCh38))=rep(FALSE, length(seqinfo(test_GRanges_GRCh38)))
seqlengths=seqlengths(seqinfo(gtf_GRCh38))[as.character(7)]
names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(test_GRanges_GRCh38)=seqlengths
genome(test_GRanges_GRCh38) <- "hg38"

Vista_GRanges_GRCh38<- test_GRanges_GRCh38 #chr7 156791088-156791875

#Convert to mm39 

ChainFilePath='/projappl/project_2006203/liftOver/hg38ToMm39.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
Vista_GRanges_Mm39=unlist(rtracklayer::liftOver(Vista_GRanges_GRCh38, ch))

Vista_GRanges_Mm39_single_range=reduce(Vista_GRanges_Mm39, drop.empty.ranges = TRUE, min.gapwidt=10)
isCircular(seqinfo(Vista_GRanges_Mm39_single_range))=rep(FALSE, length(seqinfo(Vista_GRanges_Mm39_single_range)))
seqlengths=seqlengths(seqinfo(bsg))[unique(as.character(seqnames(Vista_GRanges_Mm39)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(Vista_GRanges_Mm39_single_range)=seqlengths
genome(Vista_GRanges_Mm39_single_range) <- "mm39"

#Convert to mm10

ChainFilePath='/projappl/project_2006203/liftOver/mm39ToMm10.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
Vista_GRanges_mm10=unlist(rtracklayer::liftOver(Vista_GRanges_Mm39_single_range, ch))


isCircular(seqinfo(Vista_GRanges_mm10))=rep(FALSE, length(seqinfo(Vista_GRanges_mm10)))
seqlengths=seqlengths(seqinfo(bsg_mm10))[unique(as.character(seqnames(Vista_GRanges_mm10)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(Vista_GRanges_mm10)=seqlengths
genome(Vista_GRanges_mm10) <- "mm10"



export.bed(Vista_GRanges_Mm39, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/VISTA_SHH_enhancer_hs2496_align_to_human.bed"))
export.bed(Vista_GRanges_Mm39_single_range, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/VISTA_SHH_enhancer_hs2496.bed"))

#Weak ETS sites from literature

ETS_sites_path="/projappl/project_2006203/TFBS/ZRS_enhancer/ETS_sites/mm39/"
ETS_sites=dir(ETS_sites_path)
ETS_matches=ETS_sites[grep("match", ETS_sites)]
ETS_sites=ETS_sites[-grep("match", ETS_sites)]

ETS_sites_GRanges=GRangesList()
for(es in ETS_sites){

  ETS_sites_GRanges[[gsub(".bed", "",es)]]=import.bed(paste0(ETS_sites_path, es))
    
}

ETS_matches_GRanges=GRangesList()
for(es in ETS_matches){
  
  ETS_matches_GRanges[[gsub("_match.bed", "",es)]]=import.bed(paste0(ETS_sites_path, es))
  
}

ETS_sites_GRanges_final=unlist(ETS_sites_GRanges[unique(findOverlaps( unlist(ETS_matches_GRanges), unlist(ETS_sites_GRanges), type="within")@to)][names(ETS_sites_GRanges)])
# ETS-1     chr5 29520450-29520467      - |        <NA>      1000
# ETS-3     chr5 29520068-29520084      - |        <NA>      1000
# ETS-2     chr5 29520174-29520190      - |        <NA>      1000
# ETS-5     chr5 29519963-29519979      - |        <NA>      1000
# ETS-4     chr5 29520007-29520023      - |        <NA>      1000
# ETS-A_from_Lim_etal_2022     chr5 29520315-29520348      - |        <NA>      1000
# ETV-A     chr5 29520329-29520345      - |        <NA>      1000
ETS_matches_GRanges_final=unlist(ETS_matches_GRanges)[unique(findOverlaps( unlist(ETS_matches_GRanges), unlist(ETS_sites_GRanges), type="within")@from)][unique(findOverlaps( unlist(ETS_matches_GRanges), unlist(ETS_sites_GRanges), type="within")@to)]
# ETS-1     chr5 29520456-29520462      - |        <NA>      1000
# ETS-2_ETS-3     chr5 29520073-29520079      - |        <NA>      1000 ETS3
# ETS-2_ETS-3     chr5 29520179-29520185      - |        <NA>      1000 ETS2
# ETS-4_ETS-5     chr5 29519968-29519974      - |        <NA>      1000 ETS5
# ETS-4_ETS-5     chr5 29520012-29520018      - |        <NA>      1000 ETS4
# ETS-A_from_Lim_etal_2022     chr5 29520328-29520335      - |        <NA>      1000
# ETV-A     chr5 29520334-29520340      - |        <NA>      1000
names(ETS_matches_GRanges_final)=names(ETS_sites_GRanges_final)

isCircular(seqinfo(ETS_matches_GRanges_final))=rep(FALSE, length(seqinfo(ETS_matches_GRanges_final)))
seqlengths=seqlengths(seqinfo(bsg))[unique(as.character(seqnames(ETS_matches_GRanges_final)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(ETS_matches_GRanges_final)=seqlengths
genome(ETS_matches_GRanges_final) <- "mm39"

isCircular(seqinfo(ETS_sites_GRanges_final))=rep(FALSE, length(seqinfo(ETS_sites_GRanges_final)))
seqlengths=seqlengths(seqinfo(bsg))[unique(as.character(seqnames(ETS_sites_GRanges_final)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(ETS_sites_GRanges_final)=seqlengths
genome(ETS_sites_GRanges_final) <- "mm39"

#findOverlaps(ETS_matches_GRanges_final, ETS_sites_GRanges_final)


# mouse genomes

GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/representative_motif_matches_mouse_final.Rds")) #

#Only ETS, all or only representatives

ETS=representatives[grep("*Ets*", representatives$Lambert2018_families),] #tmp
ETS_representatives=ETS[ETS$new_representative=="YES",] #168

ETS_monomers=ETS_representatives[which(!ETS_representatives$experiment=="CAP-SELEX"),]



dimer_representatives=rep_motifs[which(rep_motifs$experiment=="CAP-SELEX"),] #596

#genome(Vista_GRanges_GRCh38) <- "GRCh38"

for(i in 1:nrow(dimer_representatives) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
 
  name=dimer_representatives$ID[i]
  
  TFBS_GRanges=GR_list[[name]]
  
  print(name)
  print(length(findOverlaps( TFBS_GRanges,Vista_GRanges_Mm39_single_range, ignore.strand=TRUE)))
  
  hits=findOverlaps( TFBS_GRanges,Vista_GRanges_Mm39_single_range, ignore.strand=TRUE)
  if(length(hits)!=0){
  export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/bed_dimer_rep/",dimer_representatives$ID[i],".bed") )
  }

}

#All representative motif hits at ETS sites and matches and Collect to a table

GRanges_list_sites<-list()
GRanges_list_matches<-list()

GRanges_list_sites[["ETS-1"]]=GRanges()
GRanges_list_sites[["ETS-2"]]=GRanges()
GRanges_list_sites[["ETS-3"]]=GRanges()
GRanges_list_sites[["ETS-4"]]=GRanges()
GRanges_list_sites[["ETS-5"]]=GRanges()
GRanges_list_sites[["ETS-A_from_Lim_etal_2022"]]=GRanges()
GRanges_list_sites[["ETV-A"]]=GRanges()

GRanges_list_matches[["ETS-1"]]=GRanges()
GRanges_list_matches[["ETS-2"]]=GRanges()
GRanges_list_matches[["ETS-3"]]=GRanges()
GRanges_list_matches[["ETS-4"]]=GRanges()
GRanges_list_matches[["ETS-5"]]=GRanges()
GRanges_list_matches[["ETS-A_from_Lim_etal_2022"]]=GRanges()
GRanges_list_matches[["ETV-A"]]=GRanges()




for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
 
  name=rep_motifs$ID[i]
 
  TFBS_GRanges=GR_list[[name]]
  
  hits=findOverlaps( TFBS_GRanges,ETS_sites_GRanges_final, ignore.strand=TRUE)
  if(length(hits)!=0){
    export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/bed_rep_site_overlap/",name,".bed") )
  }
  
  hits=findOverlaps( TFBS_GRanges,ETS_matches_GRanges_final, ignore.strand=TRUE)
  if(length(hits)!=0){
    export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/bed_rep_match_overlap/",name,".bed") )
  }
  
 
  for(j in 1:length(ETS_sites_GRanges_final)){

  #length(ETS_match_GRanges_final)
  
  print(name)
  
  print(length(findOverlaps( TFBS_GRanges,ETS_sites_GRanges_final[j], ignore.strand=TRUE))) #Or just the VISTA enhancer
  
  hits=findOverlaps( TFBS_GRanges,ETS_sites_GRanges_final[j], ignore.strand=TRUE)
  if(length(hits)!=0){
    
    
    site_matches=TFBS_GRanges[hits@from]
    site_matches$name=name
    GRanges_list_sites[[j]]=c(GRanges_list_sites[[j]], site_matches)
    
    hits=findOverlaps( TFBS_GRanges,ETS_matches_GRanges_final[j], ignore.strand=TRUE)
    if(length(hits)!=0){
      match_matches=TFBS_GRanges[hits@from]
      match_matches$name=name
      GRanges_list_matches[[j]]=c(GRanges_list_matches[[j]], match_matches)
    }
    
  }
  
  
  }
  
}


#Let's sort GRanges_all by score
df_list_sites<-list()
df_list_matches<-list()

for(j in 1:length(GRanges_list_matches)){
  tmp=as.data.frame(GRanges_list_matches[[j]][order(GRanges_list_matches[[j]]$score, decreasing=TRUE)])
  tmp$family=representatives$Lambert2018_families[match(tmp$name, representatives$ID)]
  tmp$study=representatives$study[match(tmp$name, representatives$ID)]
  df_list_matches[[names(GRanges_list_matches)[j]]]=tmp
  
  tmp=as.data.frame(GRanges_list_sites[[j]][order(GRanges_list_sites[[j]]$score, decreasing=TRUE)])
  tmp$family=representatives$Lambert2018_families[match(tmp$name, representatives$ID)]
  tmp$study=representatives$study[match(tmp$name, representatives$ID)]
  df_list_sites[[names(GRanges_list_matches)[j]]]=tmp
}


save.image("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/RData/mouse_ZRS.RData")


library("openxlsx")

# Create a new workbook
wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_matches)){

  addWorksheet(wb, names(df_list_matches)[j])
  writeData(wb, names(df_list_matches)[j], x=df_list_matches[[j]], colNames=TRUE)
}

saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/SHH_enhancer_ETS_match_matches",".xlsx"), overwrite = TRUE)


wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_sites)){
  
  addWorksheet(wb, names(df_list_sites)[j])
  writeData(wb, names(df_list_sites)[j], x=df_list_sites[[j]], colNames=TRUE)
}

saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/SHH_enhancer_ETS_sites_matches",".xlsx"), overwrite = TRUE)


#Representative motif pwms 


motifs<-readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pfms.Rds")

#motifmatchr expects input PWMs to use either natural logarithm or log 2. 
#If the input is a PFM, the TFBSTools toPWM is used for making the PWM 
#with the default psueodcounts of 0.8 and base 2 logarithm. For more control of the pseudocounts, 
#simply use the toPWM function to convert your PFMs prior to calling matchMotifs.

#MOODS
#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#  --log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#  --lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

library(TFBSTools)
#example_pwms <- do.call(PWMatrixList,lapply(motifs, toPWM, 
#                                            pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))

#saveRDS(example_pwms, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")
example_pwms=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")


#Let's match motifs directly to the VISTA enhancer to find the weak monomeric ETS sites

ETS_monomer_pwms=example_pwms[ETS_monomers$ID]

#p.cutoff=5e-02 0.05
ETS_motif_ix <- matchMotifs(ETS_monomer_pwms, Vista_GRanges_Mm39_single_range, genome="mm39", p.cutoff=0.1, out="positions") #out=scores

#sort based on the max score

ETS_motif_ix_maxscore<-list()

ETS_motif_ix_lnscore<-list()

for(i in 1:length(ETS_motif_ix)){
  print(i)
  ETS_motif_ix[[i]]=ETS_motif_ix[[i]][rev(order(ETS_motif_ix[[i]]$score))]
  ETS_motif_ix_lnscore[[i]]=ETS_motif_ix[[i]]
  ETS_motif_ix_lnscore[[i]]$score=ETS_motif_ix[[i]]$score/log2(exp(1))
  ETS_motif_ix_lnscore[[i]]=ETS_motif_ix_lnscore[[i]][ETS_motif_ix_lnscore[[i]]$score>0]
  if(length(ETS_motif_ix_lnscore[[i]]) != 0){
    export.bed(ETS_motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/bed_ETS_monomer_relaxed/",names(ETS_motif_ix)[[i]],".bed") )
    export.bedGraph(ETS_motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/bed_ETS_monomer_relaxed/",names(ETS_motif_ix)[[i]],".bedGraph") )
    ETS_motif_ix_maxscore[[names(ETS_motif_ix)[[i]] ]]=as.data.frame(ETS_motif_ix_lnscore[[i]][1,])
  }
}

#Evgeny Kvon
#ZRS has functional binding sites for ETS1 (we showed that they are mutated in snakes), 
#HOXD, E-box and ETV factors. There are also SOX repressor binding sites. 

#Please find attached two mm10 lists of limb enhancers, short-range (10-200 kb range) 
#and long-range (>400 kb away). These are all confirmed limb enhancers in transgenic assays
#and were linked to genes based on multiome or capture-HiC analysis. 

#We ran some basic motif analysis based on homer and found that various 
#Hox motifs (Hoxd9, Hoxd13 etc.), Lhx motifs and other homeodomain motifs 
#are strongly enriched in distal enhancers. 
#This is still unpublished and will be part of our ZRS enhancer swap paper 
#so please also keep confidential for now. 
#It would be interesting to see if your better PWM models could detect more matches. 

# Does these contain the corresponding sequence of this Human element [ hs2496 ]

short_dist_limb_enh=read_delim("/projappl/project_2006203/TFBS/EvgenyKvon_limbEnhancers_mm10/shortdistance_10-200kb_LimbEnhancers.txt", delim="\t")
long_dist_limb_enh=read_delim("/projappl/project_2006203/TFBS/EvgenyKvon_limbEnhancers_mm10/longdist_400-2MbLimbEnhancers.txt", delim="\t")
long_dist_limb_enh[which(long_dist_limb_enh$VistaEnh=="hs2496"),] #mm39 chr5:29519879-29520665 

#mm10 chr5:29314274-29316273 #Vista_GRanges_mm10: chr5 29314882-29315667      *

GRanges_short_dist_limb_enh=GRanges(short_dist_limb_enh) #142
GRanges_long_dist_limb_enh=GRanges(long_dist_limb_enh) #26

isCircular(seqinfo(GRanges_short_dist_limb_enh))=rep(FALSE, length(seqinfo(GRanges_short_dist_limb_enh)))
seqlengths=seqlengths(seqinfo(bsg_mm10))[unique(as.character(seqnames(GRanges_short_dist_limb_enh)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(GRanges_short_dist_limb_enh)=seqlengths[names(seqlengths(GRanges_short_dist_limb_enh))]
genome(GRanges_short_dist_limb_enh) <- "mm10"


isCircular(seqinfo(GRanges_long_dist_limb_enh))=rep(FALSE, length(seqinfo(GRanges_long_dist_limb_enh)))
seqlengths=seqlengths(seqinfo(bsg_mm10))[unique(as.character(seqnames(GRanges_long_dist_limb_enh)))]
#names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(GRanges_long_dist_limb_enh)=seqlengths[names(seqlengths(GRanges_long_dist_limb_enh))]
genome(GRanges_long_dist_limb_enh) <- "mm10"



GRanges_list_short<-list()
GRanges_list_long<-list()

for(i in 1:length(GRanges_short_dist_limb_enh)){
  #print(GRanges_short_dist_limb_end[i]$VistaEnh)
  GRanges_list_short[[GRanges_short_dist_limb_enh[i]$VistaEnh]]=GRanges()
}


for(i in 1:length(GRanges_long_dist_limb_enh)){
  #print(name(se))
  GRanges_list_long[[GRanges_long_dist_limb_enh[i]$VistaEnh]]=GRanges()
  
}

#Convert to mm10

ChainFilePath='/projappl/project_2006203/liftOver/mm39ToMm10.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)





for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  
  name=rep_motifs$ID[i]
  
  TFBS_GRanges=unlist(rtracklayer::liftOver(GR_list[[name]], ch))
  
  
  isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  seqlengths=seqlengths(seqinfo(bsg_mm10))[unique(as.character(seqnames(TFBS_GRanges)))]
  seqlengths(TFBS_GRanges)=seqlengths[seqnames(seqinfo(TFBS_GRanges))]
  genome(TFBS_GRanges) <- "mm10"
  
  
  for(j in 1:length(GRanges_short_dist_limb_enh)){
    
    #length(ETS_match_GRanges_final)
    
    #print(name)
    
    #print(length(findOverlaps( TFBS_GRanges,GRanges_short_dist_limb_enh[j], ignore.strand=TRUE))) #Or just the VISTA enhancer
    
    hits=findOverlaps( TFBS_GRanges,GRanges_short_dist_limb_enh[j], ignore.strand=TRUE)
    if(length(hits)!=0){
      
      
      matches=TFBS_GRanges[hits@from]
      matches$name=name
      GRanges_list_short[[j]]=c(GRanges_list_short[[j]], matches)
      
      
    }
  }
  
  for(j in 1:length(GRanges_long_dist_limb_enh)){
    
    #length(ETS_match_GRanges_final)
    
    #print(name)
    
    #print(length(findOverlaps( TFBS_GRanges,GRanges_long_dist_limb_enh[j], ignore.strand=TRUE))) #Or just the VISTA enhancer
    
    hits=findOverlaps( TFBS_GRanges,GRanges_long_dist_limb_enh[j], ignore.strand=TRUE)
    if(length(hits)!=0){
      matches=TFBS_GRanges[hits@from]
      matches$name=name
      GRanges_list_long[[j]]=c(GRanges_list_long[[j]], matches)
    }
  }
}


#Let's sort GRanges_all by score
df_list_sites<-list()
df_list_matches<-list()

for(j in 1:length(GRanges_list_matches)){
  df_list_matches[[names(GRanges_list_matches)[j]]]=as.data.frame(GRanges_list_matches[[j]][order(GRanges_list_matches[[j]]$score, decreasing=TRUE)])
  
  df_list_sites[[names(GRanges_list_matches)[j]]]=as.data.frame(GRanges_list_sites[[j]][order(GRanges_list_sites[[j]]$score, decreasing=TRUE)])
}

library("openxlsx")

# Create a new workbook
wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_matches)){
  
  addWorksheet(wb, names(df_list_matches)[j])
  writeData(wb, names(df_list_matches)[j], x=df_list_matches[[j]], colNames=TRUE)
}

saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/SHH_enhancer_ETS_match_matches",".xlsx"), overwrite = TRUE)


wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_sites)){
  
  addWorksheet(wb, names(df_list_sites)[j])
  writeData(wb, names(df_list_sites)[j], x=df_list_sites[[j]], colNames=TRUE)
}

saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/mouse/SHH_enhancer_ETS_sites_matches",".xlsx"), overwrite = TRUE)





# motif_ix <- matchMotifs(example_pwms, Vista_GRanges_GRCh38, genome="hg38", p.cutoff=5e-02, out="positions") #out=scores
# 
# #sort based on the max score
# 
# motif_ix_maxscore<-list()
# 
# motif_ix_lnscore<-list()
# 
# for(i in 1:length(motif_ix)){
#   print(i)
#   motif_ix[[i]]=motif_ix[[i]][rev(order(motif_ix[[i]]$score))]
#   motif_ix_lnscore[[i]]=motif_ix[[i]]
#   motif_ix_lnscore[[i]]$score=motif_ix[[i]]$score/log2(exp(1))
#   motif_ix_lnscore[[i]]=motif_ix_lnscore[[i]][motif_ix_lnscore[[i]]$score>=1]
#   if(length(motif_ix_lnscore[[i]]) != 0){
#     export.bed(motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed_relaxed/",names(motif_ix)[[i]],".bed") )
#     motif_ix_maxscore[[names(motif_ix)[[i]] ]]=as.data.frame(motif_ix_lnscore[[i]][1,])
#   }
# }
# 
# best_motif_matches=do.call(rbind, motif_ix_maxscore)
# best_motif_matches$pfm=names(motif_ix_maxscore)
# 
# #sort based on motif match
# best_motif_matches=best_motif_matches[order(best_motif_matches$score, decreasing=TRUE),]


