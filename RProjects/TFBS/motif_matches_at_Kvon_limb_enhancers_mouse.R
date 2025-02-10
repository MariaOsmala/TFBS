
library("rtracklayer")
library("TFBSTools")
library("Biostrings")
library("motifmatchr")
library("GenomicRanges")
library("readr")
library("ggplot2")


setwd("/projappl/project_2006203/TFBS/RProjects/TFBS")

library(BSgenome.Mmusculus.UCSC.mm10)
bsg_mm10 <- BSgenome.Mmusculus.UCSC.mm10
seq_len_mm10=seqlengths(bsg_mm10)
names(seq_len_mm10)=gsub("chr", "", names(seq_len_mm10))
seqinfo(bsg_mm10)
seqlevelsStyle(bsg_mm10) <- "UCSC" #1,2,3,4,..,X,Y
saveRDS(seq_len_mm10, file="/scratch/project_2006472/RData/mm10_chrom_lens.Rds")


setwd


data_path="/scratch/project_2006203/TFBS/"

representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives[which(representatives$new_representative=="YES"),]
#length(unique(rep_motifs$ID))#1031


#gtf_GRCh38 <- readRDS("/scratch/project_2006203/TFBS/RProjects/TFBS/gtf.Rds")


# mouse genomes

GR_list=readRDS(file = paste0("/scratch/project_2006203/TFBS/", "/ATAC-seq-peaks/RData/representative_motif_matches_mouse_final.Rds")) #


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

GRanges_short_dist_limb_enh=GRanges(short_dist_limb_enh) #26
GRanges_long_dist_limb_enh=GRanges(long_dist_limb_enh) #142

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

rm(GR_list)

save.image(file = paste0("/scratch/project_2006203/TFBS/","/ATAC-seq-peaks/RData/Kvon_enhancers.RData"))

load("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/Kvon_enhancers.RData")

#saveRDS(GRanges_list_long, file = paste0("/scratch/project_2006203/TFBS/", 
#                                         "/ATAC-seq-peaks/RData/Kvon_longlist_enhancers.Rds")) #


#saveRDS(GRanges_list_short, file = paste0("/scratch/project_2006203/TFBS/", 
#                                         "/ATAC-seq-peaks/RData/Kvon_shortlist_enhancers.Rds")) #


stop()

#Let's sort GRanges_all by score
df_list_short<-list()
df_list_long<-list()

for(j in 1:length(GRanges_list_short)){
  df_list_short[[names(GRanges_list_short)[j]]]=as.data.frame(GRanges_list_short[[j]][order(GRanges_list_short[[j]]$score, decreasing=TRUE)])
  
  df_list_short[[names(GRanges_list_short)[j]]]$family=representatives$Lambert2018_families[match(df_list_short[[names(GRanges_list_short)[j]]]$name, representatives$ID)]
  
  
}

for(j in 1:length(GRanges_list_long)){
  df_list_long[[names(GRanges_list_long)[j]]]=as.data.frame(GRanges_list_long[[j]][order(GRanges_list_long[[j]]$score, decreasing=TRUE)])
  
  df_list_long[[names(GRanges_list_long)[j]]]$family=representatives$Lambert2018_families[match(df_list_long[[names(GRanges_list_long)[j]]]$name, representatives$ID)]
}





library("openxlsx")

# Create a new workbook
wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_short)){
  
  addWorksheet(wb, names(df_list_short)[j] )
  writeData(wb, names(df_list_short)[j], x=df_list_short[[j]], colNames=TRUE)
}

saveWorkbook(wb,"/projappl/project_2006203/TFBS/EvgenyKvon_limbEnhancers_mm10/shortdistance_10-200kb_LimbEnhancers.xlsx", overwrite = TRUE)


wb = createWorkbook()

# Add a worksheet
for(j in 1:length(df_list_long)){
  
  addWorksheet(wb, names(df_list_long)[j])
  writeData(wb, names(df_list_long)[j], x=df_list_long[[j]], colNames=TRUE)
}


saveWorkbook(wb,"/projappl/project_2006203/TFBS/EvgenyKvon_limbEnhancers_mm10/longdist_400-2MbLimbEnhancers.xlsx", overwrite = TRUE)


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


