#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))


library("rtracklayer")
library("TFBSTools")
library("Biostrings")
library("motifmatchr")
library("GenomicRanges")
#library("BSgenome.Hsapiens.UCSC.hg38")
#library("BSgenome.Hsapiens.NCBI.GRCh38")


representatives=read.table("/projappl/project_2006203/TFBS/PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives[which(representatives$new_representative=="YES"),]

gtf_GRCh38 <- readRDS("/scratch/project_2006203/TFBS/RProjects/TFBS/gtf.Rds")

#REX element
#rex_hg19=BiocIO::import("/scratch/project_2006203/TFBS/Results/REX/REX_hg19.bed") why does not work


#This is hg19, convert to GRCh38
rex_hg19=GRanges(seqnames = Rle( "chr16", rep(1, 1) ),
                     ranges = IRanges(start=51658717, end = 51659483) )

genome(rex_hg19) <- "GRCh37"

ChainFilePath='/projappl/project_2006203/liftOver/hg19ToHg38.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
rex_GRCh38=unlist(rtracklayer::liftOver(rex_hg19, ch))
isCircular(seqinfo(rex_GRCh38))=rep(FALSE, length(seqinfo(rex_GRCh38)))
seqlengths=seqlengths(seqinfo(gtf_GRCh38))[as.character(16)]
names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(rex_GRCh38)=seqlengths
genome(rex_GRCh38) <- "GRCh38"

export.bed(rex_GRCh38, paste0("/scratch/project_2006203/TFBS/Results/REX/REX_GRCh38.bed"))


#hs72 chr16 51657810 51658716
hs72_hg19=GRanges(seqnames = Rle( "chr16", rep(1, 1) ),
                 ranges = IRanges(start=51657810, end = 51658716) )

genome(hs72_hg19) <- "GRCh37"
hs72_GRCh38=unlist(rtracklayer::liftOver(hs72_hg19, ch))
isCircular(seqinfo(hs72_GRCh38))=rep(FALSE, length(seqinfo(hs72_GRCh38)))
seqlengths=seqlengths(seqinfo(gtf_GRCh38))[as.character(16)]
names(seqlengths)=paste0("chr", names(seqlengths))
seqlengths(hs72_GRCh38)=seqlengths
genome(hs72_GRCh38) <- "GRCh38"

export.bed(hs72_GRCh38, paste0("/scratch/project_2006203/TFBS/Results/REX/hs72_GRCh38.bed"))


region=resize(reduce(c(hs72_GRCh38,rex_GRCh38)), 5000, fix="center")


#PhyloP scores
phylop.gr=import("/scratch/project_2007567/phyloP/241-mammalian-2020v2.bigWig", which=region,as = "GRanges")
BiocIO::export(phylop.gr, con="/scratch/project_2006203/TFBS/Results/REX/phyloP.wig", format="wig")

significantly_conserved <- readRDS("/projappl/project_2007567/RProject/RData/significantly_conserved.Rds")

hits=findOverlaps( significantly_conserved,region, ignore.strand=TRUE)
if(length(hits)!=0){
  export.bed( significantly_conserved[hits@from], paste0("/scratch/project_2006203/TFBS/Results/REX/significant_conserved.bed") )
}

#Conserved adjacent part?
#chr16:51624738-51624949

#GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds")) #
GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")) #


all_sites=GRanges()
for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  name=rep_motifs$ID[i]
  
  TFBS_GRanges=GR_list[[name]]
  
  hits=findOverlaps( TFBS_GRanges,rex_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
    print(i)
    export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/REX/bed_rep_site_overlap/",name,".bed") )
 
    site_matches=TFBS_GRanges[hits@from]
    df=mcols(site_matches)
    df$name=name
    mcols(site_matches)=df
    all_sites=c(all_sites, site_matches)

     

  }
    
}



df_list_sites<-list()

tmp=as.data.frame(all_sites[order(all_sites$score, decreasing=TRUE)])
tmp$family=representatives$Lambert2018_families[match(tmp$name, representatives$ID)]
tmp$study=representatives$study[match(tmp$name, representatives$ID)]

save.image("/scratch/project_2006203/TFBS/Results/REX/RData/human_REX.RData")


library("openxlsx")

# Create a new workbook
wb = createWorkbook()

# Add a worksheet

addWorksheet(wb, "REX_human_GRCh38")
writeData(wb, "REX_human_GRCh38", x=tmp, colNames=TRUE)


saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/REX/REX_matches_human_GRCh38",".xlsx"), overwrite = TRUE)


wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_sites)){
  
  addWorksheet(wb, names(df_list_sites)[j])
  writeData(wb, names(df_list_sites)[j], x=df_list_sites[[j]], colNames=TRUE)
}

saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/SHH_enhancer_ETS_sites_matches",".xlsx"), overwrite = TRUE)



#human and aligned mouse genomes




#Only ETS

ETS=representatives[grep("*Ets*", representatives$Lambert2018_families),] #tmp
ETS_representatives=ETS[ETS$new_representative=="YES",] #168

ETS_monomers=ETS_representatives[which(!ETS_representatives$experiment=="CAP-SELEX"),]


#What represents ETS1
#represented_by_representatives <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/represented_by_representatives.RDS")
#represented_by_representatives <- readRDS("/scratch/project_2006203/TFBS/RProjects/TFBS/RData/represented_by_representatives.RDS")

#indices <- which(sapply(represented_by_representatives, function(x) any(grepl("ETS1", x))))
#Answer: ELK1_HT-SELEX_TCGGAA20NAGT_AG_ACCGGAAGTN_1_2

#which represents the representative

# 
# ETS_nonrepresentatives=ETS[ETS$new_representative=="NO",] #334
# NOTHING ADDED TO THE REPRESENTATIVE SET
# for(i in 1:nrow(ETS_nonrepresentatives)){
#   indices <- which(sapply(represented_by_representatives, function(x) any(grepl(ETS_nonrepresentatives$ID[1], x))))
#   
#   print(representatives$Lambert2018_families[which(representatives$ID %in% names(indices))])
#   
#   #If none of the representatives are not in the representative set already defined, add one to the representative set
#   if( length(which(names(indices) %in% ETS_representatives$ID==TRUE))==0){
#     ETS_representatives=cbind(ETS_representatives, representatives[which(representatives$ID %in% names(indices)[1])])
#   }
# }

for(i in 1:nrow(ETS_representatives) ){
 
  print(i)
  
  #i=1
  
  name=ETS_representatives$ID[i]
  
  TFBS_GRanges=GR_list[[name]]
  
  #isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  
  #seqlengths=seqlengths(seqinfo(gtf_GRCh38))[c(as.character(seq(1,22,1)), "X", "Y")]
  #names(seqlengths)=paste0("chr", names(seqlengths))
  
  #seqlevels(TFBS_GRanges)=paste0("chr",seqlevels(seqinfo(keepSeqlevels(gtf_GRCh38, c(as.character(seq(1,22,1)), "X", "Y"), pruning.mode="coarse"))))
  
  #seqlengths(TFBS_GRanges)=seqlengths
  
  #genome(TFBS_GRanges) <- "GRCh38"
  
  #print(name)
  print(length(findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE)))
  
  hits=findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
    print(name)
    export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/ETS_bed/",ETS_representatives$ID[i],".bed") )
  }
  
}




dimer_representatives=rep_motifs[which(rep_motifs$experiment=="CAP-SELEX"),] #596

genome(Vista_GRanges_GRCh38) <- "GRCh38"

for(i in 1:nrow(dimer_representatives) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  #i=1
  #threshold=rep_motifs$MOODS_threshold[i]
  name=dimer_representatives$ID[i]
  #rep_motifs$match_numbers[i]
  #if(threshold==5){
  #  TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/",name,"_top.bed"))
  #}else{
  #  TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_",threshold,"_processed/MOODS_bigbed/",name,"_top.bed"))
  #}
  TFBS_GRanges=GR_list[[name]]
  
  #isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  
  #seqlengths=seqlengths(seqinfo(gtf_GRCh38))[c(as.character(seq(1,22,1)), "X", "Y")]
  #names(seqlengths)=paste0("chr", names(seqlengths))
  
  #seqlevels(TFBS_GRanges)=paste0("chr",seqlevels(seqinfo(keepSeqlevels(gtf_GRCh38, c(as.character(seq(1,22,1)), "X", "Y"), pruning.mode="coarse"))))
  
  #seqlengths(TFBS_GRanges)=seqlengths
  
  #genome(TFBS_GRanges) <- "GRCh38"
  
  print(name)
  print(length(findOverlaps( TFBS_GRanges,Vista_GRanges_GRCh38, ignore.strand=TRUE)))
  
  hits=findOverlaps( TFBS_GRanges,Vista_GRanges_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
  export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed_dimer_rep/",dimer_representatives$ID[i],".bed") )
  }

}

#All representative motif hits and Collect to a table

GRanges_all<-GRanges()

for(i in 1:nrow(rep_motifs) ){
  #TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  
  #i=1
  #threshold=rep_motifs$MOODS_threshold[i]
  name=rep_motifs$ID[i]
  #rep_motifs$match_numbers[i]
  # if(threshold==5){
  #   TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed/",name,"_top.bed"))
  # }else{
  #   TFBS_GRanges=import.bed(paste0("/scratch/project_2006203/TFBS/Results/MOODS_Teemu_",threshold,"_processed/MOODS_bigbed/",name,"_top.bed"))
  # }
  
  TFBS_GRanges=GR_list[[name]]
  
  # isCircular(seqinfo(TFBS_GRanges))=rep(FALSE, length(seqinfo(TFBS_GRanges)))
  # 
  # seqlengths=seqlengths(seqinfo(gtf_GRCh38))[c(as.character(seq(1,22,1)), "X", "Y")]
  # names(seqlengths)=paste0("chr", names(seqlengths))
  # 
  # seqlevels(TFBS_GRanges)=paste0("chr",seqlevels(seqinfo(keepSeqlevels(gtf_GRCh38, c(as.character(seq(1,22,1)), "X", "Y"), pruning.mode="coarse"))))
  # 
  # seqlengths(TFBS_GRanges)=seqlengths
  # 
  # genome(TFBS_GRanges) <- "GRCh38"
  
  print(name)
  print(length(findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE))) #Or just the VISTA enhancer
  
  hits=findOverlaps( TFBS_GRanges,test_GRanges_GRCh38, ignore.strand=TRUE)
  if(length(hits)!=0){
    
    export.bed( TFBS_GRanges[hits@from], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed/",rep_motifs$ID[i],".bed") )
    
    tmp=TFBS_GRanges[hits@from]
    tmp$name=name
    GRanges_all=c(GRanges_all, tmp)
    
  }
  
}


#Let's sort GRanges_all by score

tmp=as.data.frame(GRanges_all[order(GRanges_all$score, decreasing=TRUE)])
write.table(tmp, quote=FALSE,file="/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/SHH_enhancer_ZRS_matches.tsv",sep="\t", row.names = FALSE)

#Representative motif pwms 

motifs<-PFMatrixList()

for(i in 1:nrow(rep_motifs)){
  print(i)
  profileMatrix=as.matrix(read.table(file=paste0("/projappl/project_2006203/TFBS/",gsub("../../", "",rep_motifs$filename[i]))))
  dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
  motifs[[rep_motifs$ID[i]]] <- PFMatrix(ID=rep_motifs$ID[i], name=rep_motifs$symbol[i], 
                                         #matrixClass="Zipper-Type", 
                                         strand="+", 
                                         bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                         tags=list(family=rep_motifs$Lambert2018.families[i], 
                                                   #species="10090", 
                                                   #tax_group="vertebrates", 
                                                   #medline="7592839", 
                                                   #type="SELEX", 
                                                   #ACC="P53762", 
                                                   #pazar_tf_id="TF0000003",
                                                   T#FBSshape_ID="11", 
                                                   #TFencyclopedia_ID="580"
                                         ),
                                         profileMatrix=profileMatrix
  )
}

#saveRDS(motifs, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pfms.Rds")
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
example_pwms <- do.call(PWMatrixList,lapply(motifs, toPWM, 
                                            pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))

#saveRDS(example_pwms, file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")
example_pwms=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/representative_pwms.Rds")


#Let's match motifs directly to the VISTA enhancer to find the weak monomeric ETS sites

ETS_monomer_pwms=example_pwms[ETS_monomers$ID]

#p.cutoff=5e-02 0.05
ETS_motif_ix <- matchMotifs(ETS_monomer_pwms, Vista_GRanges_GRCh38, genome="hg38", p.cutoff=0.1, out="positions") #out=scores

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
    export.bed(ETS_motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed_ETS_monomer_relaxed/",names(ETS_motif_ix)[[i]],".bed") )
    export.bedGraph(ETS_motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed_ETS_monomer_relaxed/",names(ETS_motif_ix)[[i]],".bedGraph") )
    ETS_motif_ix_maxscore[[names(ETS_motif_ix)[[i]] ]]=as.data.frame(ETS_motif_ix_lnscore[[i]][1,])
  }
}



motif_ix <- matchMotifs(example_pwms, Vista_GRanges_GRCh38, genome="hg38", p.cutoff=5e-02, out="positions") #out=scores

#sort based on the max score

motif_ix_maxscore<-list()

motif_ix_lnscore<-list()

for(i in 1:length(motif_ix)){
  print(i)
  motif_ix[[i]]=motif_ix[[i]][rev(order(motif_ix[[i]]$score))]
  motif_ix_lnscore[[i]]=motif_ix[[i]]
  motif_ix_lnscore[[i]]$score=motif_ix[[i]]$score/log2(exp(1))
  motif_ix_lnscore[[i]]=motif_ix_lnscore[[i]][motif_ix_lnscore[[i]]$score>=1]
  if(length(motif_ix_lnscore[[i]]) != 0){
    export.bed(motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/300000matches_SHH_enhancer_ZRS_wider_final/bed_relaxed/",names(motif_ix)[[i]],".bed") )
    motif_ix_maxscore[[names(motif_ix)[[i]] ]]=as.data.frame(motif_ix_lnscore[[i]][1,])
  }
}

best_motif_matches=do.call(rbind, motif_ix_maxscore)
best_motif_matches$pfm=names(motif_ix_maxscore)

#sort based on motif match
best_motif_matches=best_motif_matches[order(best_motif_matches$score, decreasing=TRUE),]


