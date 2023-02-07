library(rtracklayer)
library(GenomicRanges)



#Read motif matches
#dir("../../Results/MOODS_RDS/")

#files=list.files("../../Results/MOODS_RDS/",pattern = "_smaller") #399
#files=list.files("../../Results/MOODS_Vierstra_processed/MOODS_Vierstra_RDS/") #399

files=list.files("../../Results/MOODS_Teemu_processed/MOODS_RDS/") #399
files=files[grep("top", files)]

GR_list=GRangesList()

for(file in files){
  print(file)
  #GR_list=c(GR_list, readRDS(paste0("../../Results/MOODS_RDS/",file)) )
  #GR_list=c(GR_list, readRDS(paste0("../../Results/MOODS_Vierstra_processed/MOODS_Vierstra_RDS/",file)) )
  GR_list=c(GR_list, readRDS(paste0("../../Results/MOODS_Teemu_processed/MOODS_RDS/",file)) )
  
}

#saveRDS(GR_list,  file = paste0( "../RData/all_motif_matches_sbatch.Rds")) #
#saveRDS(GR_list,  file = paste0( "../RData/all_motif_matches_Vierstra.Rds")) #
saveRDS(GR_list,  file = paste0( "../RData/top_motif_matches_Teemu.Rds")) #

stop()

#The distribution of motif hits number for representatives
all_motif_matches_Vierstra <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/all_motif_matches_Vierstra.Rds")

match_numbers=lapply(all_motif_matches_Vierstra, length)

saveRDS(unlist(match_numbers),"/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Vierstra.Rds")
#match_numbers_Vierstra <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/match_numbers_Vierstra.Rds")

library(ggplot2)
# Basic histogram

df=data.frame(match_nro=unlist(match_numbers_Vierstra))

max(df$match_nro) #4063189
min(df$match_nro) #15846

ggplot(df, aes(x=match_nro)) + geom_histogram()
# Change the width of bins
pdf(file = "../Figures/Vierstra_motif-hit-nros.pdf" )
    #width, height, onefile, family, title, fonts, version,
    #paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
    #useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=match_nro)) + geom_histogram(binwidth=10000)
dev.off()
hist(unlist(lapply(all_motif_matches_Vierstra, length)))

#The distribution of motif hits number for representatives

#Which are representative motifs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))