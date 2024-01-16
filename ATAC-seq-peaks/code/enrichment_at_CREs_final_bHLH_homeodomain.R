library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
#.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))
.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new")
#library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")

#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"
#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

bHLHhomeo_all=representatives[which(representatives$Lambert2018_families %in% c("bHLH_Homeodomain","Homeodomain_bHLH")),] #38
rep_motifs=bHLHhomeo_all$ID


#grep("CTCF",rep_motifs) CTCF is in representative motif

length(unique(rep_motifs))
#38


representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/top_motif_matches_human_final.Rds")


motifs=names(representative_motif_matches)


cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
length(cCREs_list)

p_matrix=matrix(nrow=length(rep_motifs), ncol=length(cCREs_list),
                dimnames = list( rep_motifs,
                                names(cCREs_list) )) #3294 x 111

e_matrix=matrix(nrow=length(rep_motifs), ncol=length(cCREs_list),
                dimnames = list(rep_motifs,
                                names(cCREs_list) )) #3294 x 111



all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs)


N=length(all_cCREs) #all possible cCREs #435142

#m: a list of all cCREs that have a motif match#

for(hits in rep_motifs){
  print(hits)
  #hits=motifs[1]

 
  #m: the size of all cCREs that have a motif match#
  fol=findOverlaps(representative_motif_matches[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  
  m=length(unique(fol@to)) #This is now number of all cCREs that have at least one motif hit #9857
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #findOverlaps(query, subject)
    fol=findOverlaps(representative_motif_matches[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
    #col=countOverlaps(hits, cCREs_list[[cell_type]], ignore.strand=TRUE)
    
    # number of cCREs that have motif match
    
    k=length(unique(fol@to)) #This is now number of cCREs that have at least one motif hit #178
    #k=length(unique(fol@from)) This is now number of hits that overlap with cCREs (there can be multiple motifs hits in a single CRE)
    
    #m	the number of white balls in the urn.
    # number of motif matches, 300 000
    
    # N-m
    # the number of black balls in the urn.
    # the number of all possible cCREs - the number of motif matches
    
    #Nm=N-m #list of all CREs that do not have a motif This is computed in FishersExact
    
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

saveRDS(p_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_bHLH_homeo.Rds")) #
saveRDS(e_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_bHLH_homeo.Rds")) #





library("rtracklayer")
library("dbplyr")
library("dplyr")
library("tidyverse")
library("pheatmap")
library("ComplexHeatmap")
library("cluster")
library("readr")
library("grid")
library("RColorBrewer")
#library("heatmaply")
library("ggrepel") #
library("circlize")
library(gridExtra)
library(grid) # for textGrob
#install.packages("latex2exp")
#library("latex2exp")


representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

bHLHhomeo_all=representatives[which(representatives$Lambert2018_families %in% c("bHLH_Homeodomain","Homeodomain_bHLH")),] #38
rep_motifs=bHLHhomeo_all$ID

data_path="/scratch/project_2006203/TFBS/"
p_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs.Rds")) #
e_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs.Rds")) #

p_matrix=p_matrix[which(rownames(p_matrix) %in% rep_motifs),]
e_matrix=e_matrix[which(rownames(e_matrix) %in% rep_motifs),]

#source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
source("../code/heatmap_motor.R")

#Cell type names mapping

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

Adult_Celltypes_mapping <- read_delim("../CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)


#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"

#Cell type groups from Zhang et al. Figure 1

cell_type_groups<-list()
cell_type_groups[["Epithelial 1"]]=c("Parietal Cell",
                                     "Chief Cell",
                                     "Pancreatic Acinar Cell",
                                     "Ductal Cell (Pancreatic)",
                                     "Gastric Neuroendocrine Cell",
                                     "Foveolar Cell",
                                     "Hepatocyte")
cell_type_groups[["Stromal"]]=c("Mesothelial Cell", #Divide this into two groups
                                "Adipocyte",
                                "Luteal Cell (Ovarian)",
                                "Schwann Cell (General)",
                                "Pericyte (General) 1",
                                "Fibroblast (Peripheral Nerve)",
                                "Smooth Muscle (Uterine)",
                                "Smooth Muscle (General)",
                                "Vascular Smooth Muscle 2",
                                "Vascular Smooth Muscle 1",
                                "Pericyte (General) 2",
                                "Melanocyte",
                                "Smooth Muscle (Esophageal Muscularis) 3", #This would be the first of the second group
                                "Fibroblast (Gastrointestinal)",
                                "Cardiac Fibroblasts",
                                "Fibroblast (General)",
                                "Fibroblast (Epithelial)",
                                "Fibroblast (Sk Muscle Associated)",
                                "Peripheral Nerve Stromal",
                                "Satellite Cell")

cell_type_groups[["Smooth Muscle"]]=c("Smooth Muscle (Esophageal Muscularis) 2",
                                      "Smooth Muscle (Esophageal Muscularis) 1",
                                      "Smooth Muscle (GE Junction)",
                                      "Smooth Muscle (Colon) 2",
                                      "Smooth Muscle (Colon) 1",
                                      "Smooth Muscle (General Gastrointestinal)",
                                      "Smooth Muscle (Esophageal Mucosal)",
                                      "Smooth Muscle (Vaginal)")

cell_type_groups[["Myocyte"]]=c("Atrial Cardiomyocyte",
                                "Ventricular Cardiomyocyte",
                                "Type II Skeletal Myocyte",
                                "Type I Skeletal Myocyte")

cell_type_groups[["Mural"]]=c("Pericyte (General) 3",
                              "Cardiac Pericyte 1",
                              "Pericyte (Esophageal Muscularis)",
                              "Cardiac Pericyte 2",
                              "Pericyte (General) 4",
                              "Cardiac Pericyte 3",
                              "Fibroblast (Liver Adrenal)",
                              "Cardiac Pericyte 4"
)

cell_type_groups[["Adrenal Cortical"]]=c("Zona Fasciculata Cortical Cell",
                                         "Transitional Zone Cortical Cell",
                                         "Zona Glomerulosa Cortical Cell",
                                         "Cortical Epithelial-like"
)

cell_type_groups[["Immune Lymp"]]=c("Mast Cell",
                                    "Naive T cell",
                                    "T lymphocyte 2 (CD4+)",
                                    "T Lymphocyte 1 (CD8+)",
                                    "Natural Killer T Cell",
                                    "Memory B Cell",
                                    "Plasma Cell"
)

cell_type_groups[["Immune Myelo"]]=c("Macrophage (General,Alveolar)",
                                     "Macrophage (General)",
                                     "Microglia")

cell_type_groups[["Endothelial"]]=c("Endothelial Cell (Myocardial)",
                                    "Endothelial Cell (General) 1",
                                    "Endothelial Cell (General) 2",
                                    "Lymphatic Endothelial Cell",
                                    "Alveolar Capillary Endothelial Cell",
                                    "Endocardial Cell",
                                    "Endothelial (Exocrine Tissues)",
                                    "Endothelial Cell (General) 3",
                                    "Blood Brain Barrier Endothelial Cell")

cell_type_groups[["Brain glia"]]=c("Oligodendrocyte",
                                   "Oligodendrocyte Precursor",
                                   "Astrocyte 1",
                                   "CNS,Enteric Neuron",
                                   "Astrocyte 2"
)

cell_type_groups[["Brain neuron"]]=c("GABAergic Neuron 2",
                                     "GABAergic Neuron 1",
                                     "Glutamatergic Neuron 2",
                                     "Glutamatergic Neuron 1")

cell_type_groups[["Islet"]]=c("Pancreatic Alpha Cell 2",
                              "Pancreatic Beta Cell 2",
                              "Pancreatic Delta,Gamma cell",
                              "Pancreatic Beta Cell 1",
                              "Pancreatic Alpha Cell 1")



cell_type_groups[["Epithelial 2"]]=c("Alveolar Type 1 (AT1) Cell",
                                     "Alveolar Type 2 (AT2) Cell",
                                     "Alveolar Type 2,Immune",
                                     "Club Cell",
                                     "Thyroid Follicular Cell",
                                     "Cilliated Cell")

cell_type_groups[["Epithelial 3"]]=c("Keratinocyte 1",
                                     "Esophageal Epithelial Cell",
                                     "Keratinocyte 2",
                                     "Basal Epidermal (Skin)",
                                     "Airway Goblet Cell",
                                     "Eccrine Epidermal (Skin)",
                                     "Mammary Epithelial",
                                     "Basal Epithelial (Mammary)",
                                     "Granular Epidermal (Skin)",
                                     "Mammary Luminal Epithelial Cell 2",
                                     "Mammary Luminal Epithelial Cell 1",
                                     "Myoepithelial (Skin)")


cell_type_groups[["GI Epithelial"]]=c("Paneth Cell",
                                      "Enterochromaffin Cell",
                                      "Tuft Cell",
                                      "Small Intestinal Goblet Cell",
                                      "Small Intestinal Enterocyte",
                                      "Colonic Goblet Cell",
                                      "Colon Epithelial Cell 1",
                                      "Colon Epithelial Cell 2",
                                      "Colon Epithelial Cell 3")


length(cell_type_groups) #15, in the paper they mention 17

scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"
data_path="/scratch/project_2006203/"


#bHLH-homeodomain TF-pairs
#Which motifs are highlighted in the heatmap, do they include non-representative
#OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420 YES 6N
#OLIG2_NKX6-1_TGGCCG40NCTTT_YVIII_NMCATATGNNNNNNTAATTAN_m2_c3b0_short_20230420 NO 6N OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420"
#BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNTAATTAN_m1_c3b0_short_20230420 YES 4N

#OLIG2_OTP_TGGCCG40NCTTT_YPIII_NTAATTANNNNCATATGN_m1_c3b0_short_20230420 NO 4N BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNTAATTAN_m1_c3b0_short_20230420"
#PTF1A_NKX6-1_TCTGCC40NCAC_YVIII_NCASCTGNNNTMATTAN_m1_c3b0_short_20230420 YES 3N
#OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNCATATGN_m1_c3b0_short_20230420 NO 4N BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNTAATTAN_m1_c3b0_short_20230420"
#BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNNTAATTAN_m1_c3b0_short_20230420 NO 5N ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420"  NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420"

#OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0_short_20230420 NO 5N ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420"  NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420" 
#NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420 YES 5N
#ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420 YES 5N
#ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNNCAKMTGN_m1_c3b0_short_20230420 YES 6N

#NHLH1_HOXD8_TGAGAA40NGCG_YZIII_NCASCTGNNNNNNNNYAATTRN_m1_c3b0_short_20230420 YES 8N
#NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NCAKSTGNNNNNNTAATTAN_m1_c3b0_short_20230420 YES 6N
#NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NCAKSTGNNNNNTAATTAN_m1_c3b0_short_20230420 NO 5N ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420"


colnames(e_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(e_matrix), Adult_Celltypes_mapping$celltype) ]
colnames(p_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(p_matrix), Adult_Celltypes_mapping$celltype) ]

#which are enriched 

enriched=which(e_matrix>0) 

enrichments=e_matrix[enriched] #-Inf to 2.7
min(enrichments, na.rm=TRUE) #2.184901e-06
max(enrichments, na.rm=TRUE) #4.473739

#which are depleted 
depleted=which(e_matrix<0)
depletions=e_matrix[depleted]
min(depletions, na.rm=TRUE) #-3.0469
max(depletions, na.rm=TRUE) #-1.704847e-06

p_matrix_depleted=-p_matrix #1031 x 11
p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]






bHLHhomeo=c(
  "OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420",
  "OLIG2_NKX6-1_TGGCCG40NCTTT_YVIII_NMCATATGNNNNNNTAATTAN_m2_c3b0_short_20230420",
  "BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNTAATTAN_m1_c3b0_short_20230420",
  "OLIG2_OTP_TGGCCG40NCTTT_YPIII_NTAATTANNNNCATATGN_m1_c3b0_short_20230420",
  "PTF1A_NKX6-1_TCTGCC40NCAC_YVIII_NCASCTGNNNTMATTAN_m1_c3b0_short_20230420",
  "OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNCATATGN_m1_c3b0_short_20230420", 
  "BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NCATATGNNNNNTAATTAN_m1_c3b0_short_20230420", 
  "OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0_short_20230420",
  "NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420",
  "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNCAKMTGN_m1_c3b0_short_20230420",
  "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNNCAKMTGN_m1_c3b0_short_20230420",
  "NHLH1_HOXD8_TGAGAA40NGCG_YZIII_NCASCTGNNNNNNNNYAATTRN_m1_c3b0_short_20230420",
  "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NCAKSTGNNNNNNTAATTAN_m1_c3b0_short_20230420",
  "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NCAKSTGNNNNNTAATTAN_m1_c3b0_short_20230420")

bHLHhomeo_metadata=representatives[match( bHLHhomeo,representatives$ID),c("ID","Lambert2018_families", "study","new_representative", "type")] #15

final_matrix=p_matrix_depleted[ which(rownames(p_matrix_depleted) %in% bHLHhomeo), ]

#remove uninteresting cells
final_matrix=final_matrix[,-which( colSums((abs(final_matrix) < 50)) == nrow(final_matrix))]

#Convert the motif names (row labels) name + 5th column, Ns shortened
my_list=strsplit(rownames(final_matrix), "_")
# Find the length of the longest vector
max_length <- max(sapply(my_list, length))

# Pad shorter vectors with NA
padded_list <- lapply(my_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert the list to a data frame
row_name_info <- do.call(rbind, padded_list)
row_name_info[is.na(row_name_info[,9]),5]=row_name_info[is.na(row_name_info[,9]),6]
row_name_info=as.data.frame(row_name_info[,-seq(6,10,1)])
row_name_info$ID=rownames(final_matrix)

# Regular expression for a sequence of 'N's in the middle
#pattern <- "(?<=.)(N+)(?=.+$)"
#pattern <- "(?<!^)(N+)(?!$)"
pattern <- "(?<!^)(NN+)(?!$)"

# Find positions of the pattern
matches <- gregexpr(pattern, row_name_info[,5], perl = TRUE)
matches_info=as.data.frame(cbind(row_name_info[,5],do.call(rbind, lapply(matches,attributes))))

replacement=paste0(unlist(matches_info$capture.length) ,"N")

matches_info$new_name=""

for(i in 1:length(row_name_info[,5])){
  matches_info$new_name[i]=gsub(pattern, replacement[i], row_name_info[i,5], perl = TRUE)
}

row_name_info$new_name=matches_info$new_name

rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#Rotate column labels to that they can be read from up to down
#Add protein family info

# library("RColorBrewer")
# paletteLength=99 #100
# color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                           "RdBu")))(paletteLength)
# 
# color = colorRampPalette(brewer.pal(n = 7, name =
#                                           "Reds"))(paletteLength)
# 

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
# myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))
# 
# 
# 
# pheatmap(final_matrix, color=color,clustering_distance_rows = "correlation",
#          clustering_distance_cols = "correlation",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
#          drop_levels=TRUE,show_rownames=TRUE,
#          row_names_side="left", column_names_side="top",
#          show_colnames=TRUE,legend=TRUE, breaks= myBreaks,angle_col = "45", 
#          
#          treeheight_col = 0, treeheight_row = 0) #, cutree_rows=6,


fontsize=12

#
#row_name_colors <- ifelse(shorter_names_subset$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420"), "red", "black")
#Highlight representatives with a bold color
bold=representatives$new_representative[match(row_name_info$ID, representatives$ID)]
row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")

Protein_family_colors=c("#9d4977",
                        "#45bf52",
                        "#8d50ca",
                        "#74bf3e",
                        "#d975e2",
                        "#a2b937",
                        "#4e63d1",
                        "#c6b13d",
                        "#a782e8",
                        "#519531",
                        "#d247aa",
                        "#38c481",
                        "#e53580",
                        "#47cebe",
                        "#d03f36",
                        "#3bbac6",
                        "#e1733a",
                        "#5588e5",
                        "#df9a37",
                        "#a541a0",
                        "#86c170",
                        "#8c4fa1",
                        "#3b772c",
                        "#d25189",
                        "#4e995b",
                        "#ce435e",
                        "#6ebf92",
                        "#6e5397",
                        "#6a7a1d",
                        "#9a79bf",
                        "#9e8027",
                        "#4a71ac",
                        "#a64c28",
                        "#4facdc",
                        "#996531",
                        "#98a2e4",
                        "#6c6c2e",
                        "#db8bc0",
                        "#327243",
                        "#e6867a",
                        "#379577",
                        "#b86065",
                        "#206e54",
                        "#904455",
                        "#a1ac66",
                        "#d3a36c")


names(Protein_family_colors)=c("AP-2","TFAP","bHLH","Grainyhead","bZIP","EBF1",                   
                               "CENPB","DM",  "Ets; AT hook", "Ets", "E2F", "Forkhead",               
                               "GATA","GCM", "Nuclear receptor",        "HMG", "HMG/Sox",                 "Homeodomain",            
                               "CUT; Homeodomain",        "Homeodomain; Paired box", "Homeodomain; POU",        "POU", "HSF", "IRF",
                               "SMAD","Rel", "MADS box",                "Paired box" ,             "Prospero" ,               "Myb/SANT"   ,            
                               "MEIS","CSD", "NRF", "p53", "XPA", "RFX",
                               "RRM", "Runt","SAND","T-box","TEA", "C2H2 ZF",                
                               "C2H2 ZF; AT hook"  ,      "CCCH ZF"      ,           "Znf", "BED ZF")




pdf(file = paste0(scratch, "Figures/bHLH_homeo_selected.pdf"), width=10, 
    height=7)

ComplexHeatmap::Heatmap(final_matrix, 
                        col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                        #col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                        cluster_rows = TRUE,
                        cluster_columns = TRUE, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_rot=45, #315,
                        #row_names_gp = gpar(col = row_name_colors),
                        clustering_distance_rows = "pearson",
                        clustering_method_rows = "complete",
                        clustering_distance_columns = "pearson",
                        clustering_method_columns = "complete",
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        height = nrow(final_matrix)*unit(0.618,"cm"), 
                        width = ncol(final_matrix)*unit(1, "cm"), 
                        heatmap_legend_param=list(
                          #title = "-Log_10(P-val)",
                          #title = expression(-log[10](P-val)),
                          title = expression(paste(-log[10],"(P-value)")),
                          title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                          #title_position = "topleft",
                          grid_height = unit(8, "mm"),
                          grid_width = unit(8, "mm"),
                          #tick_length = unit(0.8, "mm"),
                          #border = NULL,
                          #at = object@levels,
                          #labels = at,
                          at = c(0, 100),
                          labels = c("0",">100 (Enriched)"),
                          labels_gp = gpar(fontsize = fontsize),
                          #labels_rot = 0,
                          #nrow = NULL,
                          #ncol = 1,
                          #by_row = FALSE,
                          legend_gp = gpar(fontsize = fontsize)
                          #legend_height = NULL,
                          #legend_width = NULL,
                          #legend_direction = c("vertical", "horizontal"),
                          # break_dist = NULL,
                          #graphics = NULL,
                          #param = NULL
                        )
)

dev.off()



#Try to use cosine-angle correlation as a distance measure
# row_order=c("OLIG2_OTP_NTAATTA4NCATATGN",     
#             "BHLHE22_EVX2_NCATATG4NTAATTAN",
#             "OLIG2_HOXB5_NCATATG6NYMATTAN", #YM=(C/T)(A/T)
#             "OLIG2_NKX6-1_NMCATATG6NTAATTAN",
#             "BHLHE22_EVX2_NCATATG5NTAATTAN",
#             "OLIG2_EVX2_NTAATTA4NCATATGN",
#             "NEUROG2_EVX2_NYMATTA5NCATATGN", 
#             "OLIG2_EVX2_NTAATTA5NCATATGN",
#             "ATOH1_EVX2_NTAATTA6NCAKMTGN", #KM=(G/T)(A/C)
#             "ATOH1_EVX2_NTAATTA5NCAKMTGN",   
#             "NHLH1_EVX2_NCAKSTG6NTAATTAN", #KS=(G/T)(G/C)
#             "NHLH1_EVX2_NCAKSTG5NTAATTAN",
#             "NHLH1_HOXD8_NCASCTG8NYAATTRN",
#             "PTF1A_NKX6-1_NCASCTG3NTMATTAN") 


#Use this order
row_order=c("OLIG2_EVX2_NTAATTA4NCATATGN", #Stromal 4
  "OLIG2_OTP_NTAATTA4NCATATGN",     #Stromal 4
"BHLHE22_EVX2_NCATATG4NTAATTAN", #Stromal 4
"PTF1A_NKX6-1_NCASCTG3NTMATTAN", #Pancreas 3
"ATOH1_EVX2_NTAATTA5NCAKMTGN",  #Astrocyte 5
"NHLH1_EVX2_NCAKSTG5NTAATTAN",#Astrocyte 5

"NEUROG2_EVX2_NYMATTA5NCATATGN", #Astrocyte 5
"OLIG2_EVX2_NTAATTA5NCATATGN",#Astrocyte 5
"BHLHE22_EVX2_NCATATG5NTAATTAN", #Astrocyte 5
"OLIG2_HOXB5_NCATATG6NYMATTAN", #YM=(C/T)(A/T)
"OLIG2_NKX6-1_NMCATATG6NTAATTAN",
"ATOH1_EVX2_NTAATTA6NCAKMTGN", #KM=(G/T)(A/C)
 "NHLH1_EVX2_NCAKSTG6NTAATTAN", #KS=(G/T)(G/C)
"NHLH1_HOXD8_NCASCTG8NYAATTRN") 

row_split = rep("6", length(row_order))
row_split[1:3] = "4"
row_split[4]="3"
row_split[5:9]="5"
row_split[14]="8"

row_split=factor(row_split, levels=c("4", "3", "5", "6", "8"))

final_matrix=final_matrix[row_order,]

#Add the protein family and cell_type_group info

row_name_info$new_name_full=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#same order for all, do not draw family info
#representatives$Lambert2018_families[match( row_name_info$ID[match(rownames(final_matrix),row_name_info$new_name_full)], representatives$ID)]
#result=as.data.frame(do.call(rbind, strsplit(family_pair_metadata$Lambert2018_families,"_")))
#names(result)=c("X1", "X2")

#Stromal first, then Islet, Then Brain glia, Then Brain neuron  
column_order=c("Fibroblast (Epithelial)",  "Peripheral Nerve Stromal", "Pancreatic Beta Cell 1", "Pancreatic Beta Cell 2", "Pancreatic Delta,Gamma cell",
"Astrocyte 1", "Astrocyte 2",  "Glutamatergic Neuron 1", "Glutamatergic Neuron 2",  "GABAergic Neuron 1", "GABAergic Neuron 2")     


final_matrix=final_matrix[, column_order]

#CEll type groups
ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]

#One needs to order the row_name_info also
row_name_info=row_name_info[match( row_order,row_name_info$new_name_full),]
bold=representatives$new_representative[match(row_name_info$ID, representatives$ID)]
row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")


# Generate a color palette
num_colors <- length(unique(ct_group_names))
palette_colors <- brewer.pal(num_colors, "Set2")  # Or any other palette

# Create a named vector of colors
annotation_colors <- setNames(palette_colors, unique(ct_group_names))

pdf(file = paste0(scratch, "Figures/bHLH_homeo_selected_representative_info_ordered.pdf"), width=12, 
    height=7)

ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                        col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                        #col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                        cluster_rows = FALSE,
                        cluster_columns = FALSE, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_rot=45, #315,
                        row_split=row_split,
                        cluster_row_slices=FALSE,
                        row_gap = unit(1.5, "mm"),
                        #row_title="spacing",
                        row_title_rot=0,
                        row_title_gp=gpar(fontface="bold"),
                        row_names_gp = gpar(fontface = row_name_fontfaces),
                        clustering_distance_rows = "pearson",
                        clustering_method_rows = "complete",
                        clustering_distance_columns = "pearson",
                        clustering_method_columns = "complete",
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        height = nrow(final_matrix)*unit(0.618,"cm"), 
                        width = ncol(final_matrix)*unit(1, "cm"), 
                        top_annotation=HeatmapAnnotation(`Cell type group`=ct_group_names,
                                                        col=list(`Cell type group`=annotation_colors),
                                                        annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                          gp=gpar(fontsize=fontsize)),
                                                        annotation_name_side="left",
                                                        #annotation_name_rot=0,
                                                        #annotation_label=c("Protein family 1","Protein family 2"),
                                                        #annotation_label="Protein\nfamily",
                                                        gp=gpar(fontsize=fontsize)),
                        heatmap_legend_param=list(
                          #title = "-Log_10(P-val)",
                          #title = expression(-log[10](P-val)),
                          title = expression(paste(-log[10],"(P-value)")),
                          title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                          #title_position = "topleft",
                          grid_height = unit(8, "mm"),
                          grid_width = unit(8, "mm"),
                          #tick_length = unit(0.8, "mm"),
                          #border = NULL,
                          #at = object@levels,
                          #labels = at,
                          at = c(0, 100),
                          labels = c("0",">100 (Enriched)"),
                          labels_gp = gpar(fontsize = fontsize),
                          #labels_rot = 0,
                          #nrow = NULL,
                          #ncol = 1,
                          #by_row = FALSE,
                          legend_gp = gpar(fontsize = fontsize)
                          #legend_height = NULL,
                          #legend_width = NULL,
                          #legend_direction = c("vertical", "horizontal"),
                          # break_dist = NULL,
                          #graphics = NULL,
                          #param = NULL
                        )
)

draw(ht_list, row_title="spacing", row_title_gp=gpar(fontface="bold"))

dev.off()

r.dend <- row_dend(ht_list)  #Extract row dendrogram
rcl.list <- row_order(ht_list)  #Extract clusters (output is a list)

rownames(final_matrix[rcl.list,])

"OLIG2_OTP_NTAATTA4NCATATGN"     
"BHLHE22_EVX2_NCATATG4NTAATTAN"
"OLIG2_HOXB5_NCATATG6NYMATTAN" #YM=(C/T)(A/T)
"OLIG2_NKX6-1_NMCATATG6NTAATTAN"
"BHLHE22_EVX2_NCATATG5NTAATTAN"
"OLIG2_EVX2_NTAATTA4NCATATGN"
"NEUROG2_EVX2_NYMATTA5NCATATGN" 
"OLIG2_EVX2_NTAATTA5NCATATGN"
"ATOH1_EVX2_NTAATTA6NCAKMTGN" #KM=(G/T)(A/C)
"ATOH1_EVX2_NTAATTA5NCAKMTGN"   
"NHLH1_EVX2_NCAKSTG6NTAATTAN" #KS=(G/T)(G/C)
"NHLH1_EVX2_NCAKSTG5NTAATTAN"
"NHLH1_HOXD8_NCASCTG8NYAATTRN"
"PTF1A_NKX6-1_NCASCTG3NTMATTAN" 



#All bHLH_homeodomain spacing motifs

#remove uninteresting cells
final_matrix=p_matrix_depleted[,-which( colSums((abs(p_matrix_depleted) < 50)) == nrow(p_matrix_depleted))]
#remove uninteresting motifs
final_matrix=final_matrix[-which( rowSums((abs(final_matrix) < 50)) == ncol(final_matrix)),] #27 x 23

#Convert the motif names (row labels) name + 5th column, Ns shortened


my_list=strsplit(rownames(final_matrix), "_")
# Find the length of the longest vector
max_length <- max(sapply(my_list, length))

# Pad shorter vectors with NA
padded_list <- lapply(my_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert the list to a data frame
row_name_info <- do.call(rbind, padded_list)
row_name_info[is.na(row_name_info[,9]),5]=row_name_info[is.na(row_name_info[,9]),6]
row_name_info=as.data.frame(row_name_info[,-seq(6,10,1)])
row_name_info$ID=rownames(final_matrix)

# Regular expression for a sequence of 'N's in the middle
#pattern <- "(?<=.)(N+)(?=.+$)"
#pattern <- "(?<!^)(N+)(?!$)"
pattern <- "(?<!^)(NN+)(?!$)"

# Find positions of the pattern
matches <- gregexpr(pattern, row_name_info[,5], perl = TRUE)
matches_info=as.data.frame(cbind(row_name_info[,5],do.call(rbind, lapply(matches,attributes))))

replacement=paste0(matches_info$capture.length, "N")

matches_info$new_name=""

for(i in 1:length(row_name_info[,5])){
  matches_info$new_name[i]=gsub(pattern, replacement[i], row_name_info[i,5], perl = TRUE)
}

#remove composites

remove_ind=which(matches_info$capture.start==-1)
matches_info=matches_info[-remove_ind,]
final_matrix=final_matrix[-remove_ind,]
row_name_info=row_name_info[-remove_ind,]
row_name_info$new_name=matches_info$new_name

rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#remove uninteresting rows and columns again
#remove uninteresting cells
final_matrix=final_matrix[,-which( colSums((abs(final_matrix) < 50)) == nrow(final_matrix))]

#Highlight representatives with a bold color
bold=representatives$new_representative[match(row_name_info$ID, representatives$ID)]

row_label_color=representatives$study[match(row_name_info$ID, representatives$ID)]
row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")
row_name_colors <- ifelse(row_label_color %in% c("fromYimeng"), "red", "black")

pdf(file = paste0(scratch, "Figures/bHLH_homeo_all.pdf"), width=12, 
    height=8)

ComplexHeatmap::Heatmap(final_matrix, 
                        col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                        #col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                        cluster_rows = TRUE,
                        cluster_columns = TRUE, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_rot=45, #315,
                        row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                        clustering_distance_rows = "pearson",
                        clustering_method_rows = "complete",
                        clustering_distance_columns = "pearson",
                        clustering_method_columns = "complete",
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        height = nrow(final_matrix)*unit(0.618,"cm"), 
                        width = ncol(final_matrix)*unit(1, "cm"), 
                        heatmap_legend_param=list(
                          #title = "-Log_10(P-val)",
                          #title = expression(-log[10](P-val)),
                          title = expression(paste(-log[10],"(P-value)")),
                          title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                          #title_position = "topleft",
                          grid_height = unit(8, "mm"),
                          grid_width = unit(8, "mm"),
                          #tick_length = unit(0.8, "mm"),
                          #border = NULL,
                          #at = object@levels,
                          #labels = at,
                          at = c(0, 100),
                          labels = c("0",">100 (Enriched)"),
                          labels_gp = gpar(fontsize = fontsize),
                          #labels_rot = 0,
                          #nrow = NULL,
                          #ncol = 1,
                          #by_row = FALSE,
                          legend_gp = gpar(fontsize = fontsize)
                          #legend_height = NULL,
                          #legend_width = NULL,
                          #legend_direction = c("vertical", "horizontal"),
                          # break_dist = NULL,
                          #graphics = NULL,
                          #param = NULL
                        )
)

dev.off()


######## Consider the logistic regression results #############

data_path="/scratch/project_2006203/TFBS/"
regression_coefs=readRDS( paste0( data_path, "ATAC-seq-peaks/RData/all_regression_coefs_maxscore.RDS"))

#Normalize the regression coefficients? Should I have normalized the coefficients already before training?
#Are the packages normalizing the scores always?

regression_coefs=scale(regression_coefs,center=TRUE, scale=TRUE)


final_matrix=regression_coefs[which(rownames(regression_coefs) %in% bHLHhomeo),]

#bHLHhomeo[which(!( bHLHhomeo %in% rownames(regression_coefs)))] #regression coefficients are only for representatives

rownames(final_matrix)=row_name_info$new_name_full[match(rownames(final_matrix),row_name_info$ID)]

#Same order of motifs and cell types as in enrichment analysis



final_matrix=final_matrix[match(row_order,rownames(final_matrix))[which(!is.na(match(row_order,rownames(final_matrix))))],]

#Consider the same cell lines as above

final_matrix=final_matrix[,which(colnames(final_matrix) %in% column_order)]
final_matrix=final_matrix[,column_order]


#final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 0.02)) ==nrow(final_matrix) )]


ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]

# Generate a color palette
num_colors <- length(unique(ct_group_names))
palette_colors <- brewer.pal(num_colors, "Set2")  # Or any other palette

# Create a named vector of colors
annotation_colors <- setNames(palette_colors, unique(ct_group_names))

pdf(file = paste0(scratch, "Figures/scaled_regcoeff_bHLH_homeo_selected_representative_ordered.pdf"), width=12, 
    height=7)

ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                                #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                                col=colorRamp2(c(-ceiling(max(abs(final_matrix))*10)/10  , 0, ceiling(max(abs(final_matrix))*10)/10 ), c("navy", "white", "firebrick3")), 
                                cluster_rows = FALSE,
                                cluster_columns = FALSE, 
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_names_side = "top",
                                row_names_side = "left",
                                column_names_rot=45, #315,
                                #row_split=row_split,
                                #cluster_row_slices=FALSE,
                                #row_gap = unit(1.5, "mm"),
                                #row_title="spacing",
                                row_title_rot=0,
                                row_title_gp=gpar(fontface="bold"),
                                #row_names_gp = gpar(fontface = row_name_fontfaces),
                                clustering_distance_rows = "pearson",
                                clustering_method_rows = "complete",
                                clustering_distance_columns = "pearson",
                                clustering_method_columns = "complete",
                                show_column_dend = FALSE,
                                show_row_dend = FALSE,
                                height = nrow(final_matrix)*unit(0.618,"cm"), 
                                width = ncol(final_matrix)*unit(1, "cm"), 
                                top_annotation=HeatmapAnnotation(`Cell type group`=ct_group_names,
                                                                 col=list(`Cell type group`=annotation_colors),
                                                                 annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                                   gp=gpar(fontsize=fontsize)),
                                                                 annotation_name_side="left",
                                                                 #annotation_name_rot=0,
                                                                 #annotation_label=c("Protein family 1","Protein family 2"),
                                                                 #annotation_label="Protein\nfamily",
                                                                 gp=gpar(fontsize=fontsize)),
                                heatmap_legend_param=list(
                                  #title = "-Log_10(P-val)",
                                  #title = expression(-log[10](P-val)),
                                  title = expression("Scaled regression coefficient"),
                                  title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                                  #title_position = "topleft",
                                  grid_height = unit(8, "mm"),
                                  grid_width = unit(8, "mm"),
                                  #tick_length = unit(0.8, "mm"),
                                  #border = NULL,
                                  #at = object@levels,
                                  #labels = at,
                                  #at = c(0, 100),
                                  #labels = c("0",">100 (Enriched)"),
                                  labels_gp = gpar(fontsize = fontsize),
                                  #labels_rot = 0,
                                  #nrow = NULL,
                                  #ncol = 1,
                                  #by_row = FALSE,
                                  legend_gp = gpar(fontsize = fontsize),
                                  #legend_height = NULL,
                                  #legend_width = NULL,
                                  legend_direction = "horizontal" #c("vertical", "horizontal"),
                                  # break_dist = NULL,
                                  #graphics = NULL,
                                  #param = NULL
                                )
)

draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()
