#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


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
library("ggplot2")
library("ggrepel") #
library("circlize")
library(gridExtra)
library(grid) # for textGrob
#install.packages("latex2exp")
#library("latex2exp")
.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new")
#install.packages("heatmaply")
library("heatmaply")
library("openxlsx")
library("gridExtra")
#source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")
source("../code/heatmap_motor.R")

#short names for TFs

#Cell type names mapping

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Adult_Celltypes_mapping <- read_delim("../CATLAS/Adult_Celltypes_mapping.csv", 
#                                       delim = ";", escape_double = FALSE, trim_ws = TRUE)


#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"



scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"
scratch="../"

data_path="/scratch/project_2006203/"

#representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames_match_numbers_MOODS_threshold.tsv", sep="\t", header=TRUE)
representatives=read.table( "../../PWMs_final/metadata_representatives_match_numbers_thresholds.tsv",sep="\t", header=TRUE) #3294 the order differs!

#representatives=read.table( "../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294 the order differs! This is the same order as in p_matrix
#representatives2=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

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

#How many bHLH_Homeodomain there are
nrow(representatives[(which(representatives$Lambert2018_families=="bHLH_Homeodomain")), c("ID", "study", "new_representative", "type")]) #25)
nrow(representatives[(which(representatives$Lambert2018_families=="Homeodomain_bHLH")),  c("ID", "study", "new_representative", "type")]) #13

bHLHhomeo_all=representatives[which(representatives$Lambert2018_families %in% c("bHLH_Homeodomain","Homeodomain_bHLH")),] #38

#which represent the non-representatives

represented_by_representatives <- readRDS("../../RProjects/TFBS/RData/represented_by_representatives.RDS")

motif=bHLHhomeo[bHLHhomeo_metadata$new_representative=="NO"][1]

for(motif in bHLHhomeo[bHLHhomeo_metadata$new_representative=="NO"]){
  #print(motif)
indices <- which(sapply(represented_by_representatives, function(x) any(grepl(motif, x))))
#can be many, Let's take the shortes, or the one with highest information content?
print(paste0(motif, " ", names(represented_by_representatives)[indices]))
}


#These are log_10(p-values)
#p_matrix=readRDS(  file = paste0(scratch, "RData/Teemu_p_matrix_human_cell_type_restricted.Rds")) #

p_matrix=readRDS(file = paste0(scratch,"RData/final_p_matrix_human_cell_type_restricted.Rds")) #

#sort the representatives table according to the rownames of p_matrix
representatives <- representatives[order(match(representatives$ID, rownames(p_matrix))), ]

#test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
test=representatives$ID
rep_motifs=test[which(representatives$new_representative=="YES")]
#rep_dimers=test[which(representatives$new_representative=="YES"& representatives$experiment=="CAP-SELEX")] #591
#rep_Yimeng=test[which(representatives$new_representative=="YES"& representatives$experiment=="CAP-SELEX"& representatives$study=="fromYimeng")] #238
length(which(rep_motifs %in% rownames(p_matrix))) #1031


#remove _pfm_composite_new and _pfm_spacing_new
#rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
#rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs)) #1031

#dim(p_matrix) 1062 x 111
#length(which(rep_dimers %in% rownames(p_matrix))) #595
#length(which(rep_Yimeng %in% rownames(p_matrix))) #238

shorter_names=representatives %>% filter(test %in% rownames(p_matrix)) %>% select(symbol, ID, experiment, 
                                              new_representative, Lambert2018_families, type)
#table(shorter_names)

shorter_names$new_motif=FALSE
shorter_names$new_motif[which(shorter_names$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420") )]=TRUE #245

#"pfm_composite_new" "pfm_spacing_new"       "20230420"

#FDR correction for the p-values?

#These are log2( (k / n) / (m / N) )
#e_matrix=readRDS(  file = paste0(scratch, "RData/Teemu_e_matrix_human_cell_type_restricted.Rds")) #
e_matrix=readRDS(  file = paste0(scratch, "RData/final_e_matrix_human_cell_type_restricted.Rds")) #

#dimnames(e_matrix)[[2]] #What column in Adult_Celltypes_mapping correspond to there
#dimnames(e_matrix)[[2]] %in% Adult_Celltypes_mapping$celltype



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

as.vector(unlist(cell_type_groups))[which(as.vector(unlist(cell_type_groups)) %in% Adult_Celltypes_mapping$`Cell type`==FALSE)]


#unlist(cell_type_groups) %in% Adult_Celltypes_mapping$`Cell type` all are contained

#Adult_Celltypes_mapping$`Cell type`[which(Adult_Celltypes_mapping$`Cell type` %in% unlist(cell_type_groups) ==FALSE)] # are all contained


#rename e_matrix and p_matrix

#v=match(dimnames(e_matrix)[[2]], Adult_Celltypes_mapping$celltype)
#test=data.frame(dimnames(e_matrix)[[2]], Adult_Celltypes_mapping$celltype[v], Adult_Celltypes_mapping$`Cell type`[v])

colnames(e_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(e_matrix), Adult_Celltypes_mapping$celltype) ]
colnames(p_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(p_matrix), Adult_Celltypes_mapping$celltype) ]

#Heatmaps showing Gene Ontology TF motifs with maximal enrichment in cell-type-restricted cCREs of selected cell types. 
#Only the most enriched TF motif in each of the previously identified motif archetypes (Vierstra et al., 2020) was selected as the representative and 
#the top 10 motifs were selected for each cell type. Color represents −log10P.

#which are enriched 

infinites=which(is.infinite(e_matrix)) #226
e_matrix[infinites]=NA
enriched=which(e_matrix>0) #1062*111=117882
#52050

#why this does not work
# r = ((enriched-1) %% ncol(e_matrix)) + 1
# c = floor((enriched-1) / ncol(e_matrix)) + 1
# e_matrix[r,c]

enrichments=e_matrix[enriched] #-Inf to 2.7
min(enrichments, na.rm=TRUE) #2.184901e-06
max(enrichments, na.rm=TRUE) #4.473739

#which are depleted 
depleted=which(e_matrix<0)
depletions=e_matrix[depleted]
min(depletions, na.rm=TRUE) #-3.0469
max(depletions, na.rm=TRUE) #-1.704847e-06

#exactly zero?
#equal=which(e_matrix==0) empty

#Convert the depleted -log_10 p-values to negative
str(depleted)

#  6 12.57143
#  6 28.57143
#  6 9.714286
# 6 5.428571
#  6 12.28571
#  6 4.857143
#  6 7.714286
# 6 3.142857
#  6  11.14286
#  6 10.28571
#  6 8
#  6 7.714286
#  6 11.14286
#  6 13.42857
#  6 8.285714


#How many of the new motifs are enriched at the cCREs specific to at least one cell-type group

table(representatives$type[representatives$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420")  ])
#20230420 pfm_composite_new   pfm_spacing_new 
#14 (spacing)              503               192 

#125 composite motifs and 120 spacing motifs were representative 

p_matrix_depleted=-p_matrix #1031 x 11
p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]

p_matrix_composites=p_matrix_depleted[representatives %>% filter( type=="pfm_composite_new" & new_representative=="YES" ) %>% pull(ID), ] #125 x 111
p_matrix_spacing=p_matrix_depleted[representatives %>% filter( type %in% c("pfm_spacing_new", "20230420") & new_representative=="YES" ) %>% pull(ID), ] #120 x 111

e_matrix_composites=e_matrix[representatives %>% filter( type=="pfm_composite_new" & new_representative=="YES" ) %>% pull(ID), ] #125 x 111
e_matrix_spacing=e_matrix[representatives %>% filter( type %in% c("pfm_spacing_new", "20230420") & new_representative=="YES" ) %>% pull(ID), ] #120 x 111


#consider also the fold change
#0.01 is the same p-value threshold as in the volcano plots
#Correct for multiple testing? Check lola
#Save these results into tables
#Add also family into

cell_type_and_group_enrichment_for_motifs <- function(boolean_matrix, cell_type_groups, representatives) {
  result=apply(boolean_matrix,1, function(x) names(which(x)) ) 
  max_length <- max(sapply(result, length))
  equalized_list <- lapply(result, function(x) {
    length(x) <- max_length  # extend the length
    x
  })
  
  cell_types_df <- do.call(rbind, equalized_list)
  
  #remove NA rows
  cell_types_df=as.data.frame(cell_types_df[-which(apply(cell_types_df, 1, function(x) length(which(is.na(x))))==ncol(cell_types_df)),])
  
  #cell-type-groups
  
  cell_types_groups_df=cell_types_df
  
  for(i in 1:nrow(cell_types_df)){
    #print(i)
    for(j in 1:ncol(cell_types_df)){
      #print(j)
      if(length(names(which(sapply(cell_type_groups, function(x) as.vector(cell_types_df[i,j]) %in% x)==TRUE)))!=0){
        cell_types_groups_df[i,j] <- names(which(sapply(cell_type_groups, function(x) as.vector(cell_types_df[i,j]) %in% x)==TRUE))
      }
    }
  }
  
  cell_types_df$family=representatives$Lambert2018_families[match(rownames(cell_types_df), representatives$ID)]
  
  cell_types_df <- cell_types_df %>% select(family, everything())
  
  cell_types_groups_df$family=representatives$Lambert2018_families[match(rownames(cell_types_groups_df), representatives$ID)]
  
  cell_types_groups_df <- cell_types_groups_df %>% select(family, everything())
  
  
  
  result_list=list()
  result_list[["cell_types"]]=cell_types_df
  result_list[["cell_types_groups"]]=cell_types_groups_df
  result_list
  
}


boolean_matrix=p_matrix_composites > -log10(0.01) & e_matrix_composites > 1
length(which(apply(boolean_matrix,1, function(x) length(which(x)) )>0)) #33



composite_strict=cell_type_and_group_enrichment_for_motifs(boolean_matrix, cell_type_groups, representatives)


boolean_matrix=p_matrix_composites > -log10(0.01) & e_matrix_composites > 0.75
length(which(apply(boolean_matrix,1, function(x) length(which(x)) )>0)) #79

composite_loose=cell_type_and_group_enrichment_for_motifs(boolean_matrix, cell_type_groups, representatives)


boolean_matrix=p_matrix_spacing > -log10(0.01) & e_matrix_spacing > 1
length(which(apply(boolean_matrix,1, function(x) length(which(x)) )>0)) #23

spacing_strict=cell_type_and_group_enrichment_for_motifs(boolean_matrix, cell_type_groups, representatives)

boolean_matrix=p_matrix_spacing > -log10(0.01) & e_matrix_spacing > 0.75
length(which(apply(boolean_matrix,1, function(x) length(which(x)) )>0)) #71

spacing_loose=cell_type_and_group_enrichment_for_motifs(boolean_matrix, cell_type_groups, representatives)

#with 33 of 125 representative composite motifs and 23 of 120 representative spacing motifs being enriched in at least one cCRE set (p-value < 0.01 and the log2 fold change > 1).
#with 79 of 125 representative composite motifs and 71 of 120 representative spacing motifs being enriched in at least one cCRE set (p-value < 0.01 and the log2 fold change > 0.75).



# Create a new workbook
wb <- createWorkbook()

# Add each data frame as a new sheet
addWorksheet(wb, "composite_strict")
writeData(wb, "composite_strict", composite_strict[["cell_types"]], rowNames = TRUE)
addWorksheet(wb, "composite_loose")
writeData(wb, "composite_loose", composite_loose[["cell_types"]], rowNames = TRUE)

addWorksheet(wb, "spacing_strict")
writeData(wb, "spacing_strict", spacing_strict[["cell_types"]], rowNames = TRUE)
addWorksheet(wb, "spacing_loose")
writeData(wb, "spacing_loose", spacing_loose[["cell_types"]], rowNames = TRUE)

# Save the workbook to a file
saveWorkbook(wb, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/enrichment_excels/cell_types.xlsx", overwrite = TRUE)

wb <- createWorkbook()

# Add each data frame as a new sheet
addWorksheet(wb, "composite_strict")
writeData(wb, "composite_strict", composite_strict[["cell_types_groups"]], rowNames = TRUE)
addWorksheet(wb, "composite_loose")
writeData(wb, "composite_loose", composite_loose[["cell_types_groups"]], rowNames = TRUE)

addWorksheet(wb, "spacing_strict")
writeData(wb, "spacing_strict", spacing_strict[["cell_types_groups"]], rowNames = TRUE)
addWorksheet(wb, "spacing_loose")
writeData(wb, "spacing_loose", spacing_loose[["cell_types_groups"]], rowNames = TRUE)

# Save the workbook to a file
saveWorkbook(wb, "/scratch/project_2006203/TFBS/ATAC-seq-peaks/enrichment_excels/cell_types_groups.xlsx", overwrite = TRUE)

#bHLH-homeodomain TF pairs where the motifs were spaced 4 bp apart 
#were heavily enriched cCRE sets from stromal(Nerve Stromal) and fibroblastic cells(Nerve Stromal) 58, 
#whereas the motifs spaced 3, 5 and 6 bp apart were enriched in cCRE 
#sets from pancreatic beta cells, spacing 3
#PTF1A_NKX6-1_TCTGCC40NCAC_YVIII_NCASCTGNNNTMATTAN_m1_c3b0_short_20230420
p_matrix_spacing["PTF1A_NKX6-1_TCTGCC40NCAC_YVIII_NCASCTGNNNTMATTAN_m1_c3b0_short_20230420", "Pancreatic Beta Cell 2"] #64.14925
e_matrix_spacing["PTF1A_NKX6-1_TCTGCC40NCAC_YVIII_NCASCTGNNNTMATTAN_m1_c3b0_short_20230420", "Pancreatic Beta Cell 2"] #0.7926459


#, astrocytes spacing 5
#NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420
p_matrix_spacing["NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420","Astrocyte 1" ] #35.62038
p_matrix_spacing["NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420","Astrocyte 2" ] #58.33645

e_matrix_spacing["NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420","Astrocyte 1" ] #0.56
e_matrix_spacing["NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0_short_20230420","Astrocyte 2" ] #0.36


#glutamatergic neurons, spacing 6 (OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420, 
#OLIG2_NKX6-1_TGGCCG40NCTTT_YVIII_NMCATATGNNNNNNTAATTAN_m2_c3b0_short_20230420 (not representative))

p_matrix_spacing["OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420", "Glutamatergic Neuron 1"] #134.2129
p_matrix_spacing["OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420", "Glutamatergic Neuron 2"] #113.8792

e_matrix_spacing["OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420", "Glutamatergic Neuron 1"] #0.617945
e_matrix_spacing["OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0_short_20230420", "Glutamatergic Neuron 2"] # 0.5856391



#rownames(p_matrix_depleted)=shorter_names$symbol
#rownames(e_matrix)=shorter_names$symbol



#pdf(file = paste0(scratch, "Figures/cell_group_heatmaps.pdf"), width = 8.27, height = 15, onefile=TRUE )

#newplot=FALSE


#The ten most enriched per each cell line
for( ct_group in names(cell_type_groups) ){
  print(ct_group)
   #ct_group="Epithelial 1"
  #Epithelial
  #Visualise the Hepatocyte, Acinar, Parietal, Chief, Foveolar, G. Neuroendo, Ductal, 
  p_matrix_depleted=-p_matrix
  p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]
   p_matrix_depleted_orig=p_matrix_depleted
   rownames(p_matrix_depleted)=shorter_names$symbol
   rownames(e_matrix)=shorter_names$symbol
   
   p_cell_group_matrix=p_matrix_depleted[, which( colnames(p_matrix_depleted) %in% cell_type_groups[[ct_group]] )]

   e_cell_group_matrix=e_matrix[, which( colnames(e_matrix) %in% cell_type_groups[[ct_group]] )]

   p_cell_group_matrix_orig=p_matrix_depleted_orig[, which( colnames(p_matrix_depleted_orig) %in% cell_type_groups[[ct_group]]  )]

   #what are the 10 most enriched motifs for these cell lines       

   most_enriched<-list()
   most_enriched_orig<-list()
   
   #Choosing the ten most enriched for each cell type maybe not sensible, as there 
   #are many more for some cell lines, try volcano plots and predictive analysis
   
   for(c in colnames(p_cell_group_matrix)){
  
     most_enriched[[c]]=names(sort(p_cell_group_matrix[,c], decreasing = TRUE)[1:10])
     most_enriched_orig[[c]]=names(sort(p_cell_group_matrix_orig[,c], decreasing = TRUE)[1:10])
   }

   
    length(unlist(most_enriched)) #70
    most_enriched=unique(unlist(most_enriched)) #37

    length(unlist(most_enriched_orig)) #70
    most_enriched_orig=unique(unlist(most_enriched_orig)) #44

    #index
    index=which(rownames(p_cell_group_matrix_orig) %in% unlist(most_enriched_orig))

    p_cell_group_matrix=p_cell_group_matrix[index,] 
    e_cell_group_matrix=e_cell_group_matrix[index,] 
    p_cell_group_matrix_orig=p_cell_group_matrix_orig[index,] 

    


#paletteLength=101
#color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                          "RdBu")))(paletteLength)

#color=colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))

fontsize=12
# Replace values greater than 100 with 100 and less than -100 with -100
p_cell_group_matrix_tmp=pmin(pmax(p_cell_group_matrix, -100), 100)

shorter_names_subset <- rownames(p_cell_group_matrix_orig) %>%
  map_df(~filter(shorter_names, ID == .x))

row_name_colors <- ifelse(shorter_names_subset$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420"), "red", "black")

#print( max(6,6*ncol(p_cell_group_matrix)/44) )
#print(2*nrow( p_cell_group_matrix)/7 )

pdf(file = paste0(scratch, "Figures/cell_group_heatmaps/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix)/44), height=4*nrow(p_cell_group_matrix)/7 )

#par(pin=c(max(6,6*ncol(p_cell_group_matrix)/44), height=2*nrow(p_cell_group_matrix)/7 ))

# if(newplot==FALSE){
#   newplot=TRUE
#   
# }else{
#   #plot.new()
# }

ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
         cluster_rows = TRUE,
         cluster_columns = TRUE, 
         show_row_names = TRUE,
         show_column_names = TRUE,
         column_names_side = "top",
         row_names_side = "left",
         row_names_gp = gpar(col = row_name_colors),
         clustering_distance_rows = "pearson",
         clustering_method_rows = "complete",
         clustering_distance_columns = "pearson",
         clustering_method_columns = "complete",
         show_column_dend = FALSE,
         show_row_dend = FALSE,
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
           at = c(-100, 0, 100),
           labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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

result=try( draw(ht, column_title=ct_group), silent=TRUE )

if( class(result)=="try-error"){
  dev.off()
  pdf(file = paste0(scratch, "Figures/cell_group_heatmaps/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix)/44), height=4*nrow(p_cell_group_matrix)/7 )
  ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                             cluster_rows = FALSE,
                             cluster_columns = FALSE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             row_names_gp = gpar(col = row_name_colors),
                             clustering_distance_rows = "pearson",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "pearson",
                             clustering_method_columns = "complete",
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
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
                               at = c(-100, 0, 100),
                               labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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
  result=try( draw(ht, column_title=ct_group), silent=TRUE )
  
}
dev.off()



#dev.off()

#Interactive heatmap


#install.packages("heatmaply")
p_cell_group_matrix_saturated=p_cell_group_matrix
saturated_threshold=100
p_cell_group_matrix_saturated[p_cell_group_matrix_saturated > saturated_threshold] <- saturated_threshold
p_cell_group_matrix_saturated[p_cell_group_matrix_saturated < -saturated_threshold] <- -saturated_threshold

#get the info for the TFs from metadata table

metadata=representatives[match(rownames(p_cell_group_matrix_orig), representatives$ID),]
#[1] "ID"                   "symbol"               "clone"                "family"               "Lambert2018.families"
  #[6] "organism"             "study"                "experiment"           "ligand"               "batch"               
  #[11] "seed"                 "multinomial"          "cycle"                "representative"       "new_representative"  
  #[16] "short"                "type"                 "comment"              "filename"             "IC"                  
  #[21] "length" "match_numbers", "MOODS_threshold"
  
  test=matrix(rep(metadata$ID,each=7), ncol=nrow(p_cell_group_matrix), byrow=FALSE)
  
  mat=t(p_cell_group_matrix_orig) #7 x 44 (63?)
  mat[]<-paste("Motif: ", matrix(rep(metadata$ID,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE) ,
               "\nFamily: ", matrix(rep(metadata$Lambert2018_families,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
               "\nStudy: ", matrix(rep(metadata$study,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
               "\n-log10 pvalue:", mat, 
               "\nlog2FC: ", t(e_cell_group_matrix),
               "\nIC: ", matrix(rep(metadata$IC,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
               "\nlength: ", matrix(rep(metadata$length,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
               "\nmatch numbers: ", matrix(rep(metadata$match_numbers,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
               "\nMOODS threshold: ", matrix(rep(metadata$threshold,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE)
  )
  
  #mat[] <- 
  #  
  #  lapply( colnames(mat), function(colname) {
  #  paste0("\nFamily ", metadata[ metadata$ID==colname, "Lambert2018.families"])} )
  # Create heatmap with heatmaply
  heatmaply(t(p_cell_group_matrix_saturated),
            col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
            legend=TRUE,
            showticklabels = T,
            ColSideColors = data.frame(Family=metadata$Lambert2018_families),
            custom_hovertext=mat,
            row_text_angle = 0,
            plot_method="plotly",
            column_text_angle = 45,
            subplot_margin = 0,
            cellnote = NULL,
            cellnote_color = "auto",
            cellnote_textposition = "middle right",
            cellnote_size = 12,
            margins = c(NA, NA, NA, NA),
            scale = NULL, # Scale the data
            symm = FALSE,  # Make the heatmap symmetric
            dendrogram = "none",  # Show dendrograms on both rows and columns
            k_row = 2,  # Show 2 clusters in row dendrogram
            k_col = 2,  # Show 2 clusters in column dendrogram
            cluster_row = FALSE,  # Cluster rows
            cluster_col = FALSE,  # Cluster columns
            xlab = "",  # X-axis label
            ylab = "",  # Y-axis label
            main = ct_group,  # Title of the plot
            colorbar = TRUE,  # Show color bar
            key.title = "-log10 p-value \n Negative (blue) indicates depletion",
            colorbar_xanchor = "left",
            colorbar_yanchor = "middle",
            colorbar_xpos = 1.025,
            colorbar_ypos = 0.25,
            #colorbar_len = 0.4,
            colorbar_thickness = 30,
            fixed=FALSE,
            fontsize_row = 14,
            fontsize_col = 14,
            side_color_colorbar_len = 0.5,
            prefix="",
            #width = NULL,
            #height = NULL,
            file=paste0("/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys/",ct_group,".html")
  )
  
  
  
  #browseURL("/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys/epithelial.html")
  #browseURL("/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys/",ct_group,".html")




}


#All new enriched, clustered heatmap
for( ct_group in names(cell_type_groups) ){
  #print(ct_group)
  #ct_group="Epithelial 1"
  #Epithelial
  #Visualise the Hepatocyte, Acinar, Parietal, Chief, Foveolar, G. Neuroendo, Ductal, 
  p_matrix_depleted=-p_matrix
  p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]
  p_matrix_depleted_orig=p_matrix_depleted
  rownames(p_matrix_depleted)=shorter_names$symbol
  rownames(e_matrix)=shorter_names$symbol
  
  p_cell_group_matrix=p_matrix_depleted[, which( colnames(p_matrix_depleted) %in% cell_type_groups[[ct_group]] )]
  
  e_cell_group_matrix=e_matrix[, which( colnames(e_matrix) %in% cell_type_groups[[ct_group]] )]
  
  p_cell_group_matrix_orig=p_matrix_depleted_orig[, which( colnames(p_matrix_depleted_orig) %in% cell_type_groups[[ct_group]]  )]
  
  
  #paletteLength=100
  #color = colorRampPalette(rev(brewer.pal(n = 7, name =
  #                                          "RdBu")))(paletteLength)
  
  #color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  #myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
  #              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))
  
  
  #pdf(file = paste0(scratch, "Figures/cell_group_whole_heatmaps/",ct_group,".pdf"), height=max(8,8*ncol(t(p_cell_group_matrix))/44), width=max(4, 4*nrow(t(p_cell_group_matrix))/7) )
  
  #pheatmap(p_cell_group_matrix, color=color, 
  #          clustering_distance_rows = "euclidean",
  #          clustering_distance_cols = "euclidean",
  #          clustering_method="complete",legend=TRUE,
  #          cluster_rows = TRUE,
  #          cluster_cols = TRUE,
  #          breaks= myBreaks, row_names_side = "left", row_dend_side = "right")
  # dev.off()
  # 
  
  
  fontsize=12
  # Replace values greater than 100 with 100 and less than -100 with -100
  p_cell_group_matrix_tmp=pmin(pmax(p_cell_group_matrix, -100), 100)
  
  shorter_names_subset <- rownames(p_cell_group_matrix_orig) %>%
    map_df(~filter(shorter_names, ID == .x))
  
  row_name_colors <- ifelse(shorter_names_subset$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420"), "red", "black")
  
  index=shorter_names_subset$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420")
  
  p_cell_group_matrix_tmp=p_cell_group_matrix_tmp[which(rownames(p_cell_group_matrix_orig) %in% shorter_names_subset$ID[ index ]),]
  
  p_cell_group_matrix_new_motifs=p_cell_group_matrix[which(rownames(p_cell_group_matrix_orig) %in% shorter_names_subset$ID[ index ]),]
  p_cell_group_matrix_orig_new_motifs=p_cell_group_matrix_orig[which(rownames(p_cell_group_matrix_orig) %in% shorter_names_subset$ID[ index ]),]
  #print( max(6,6*ncol(p_cell_group_matrix)/44) )
  #print(2*nrow( p_cell_group_matrix)/7 )
  
  #Whole heatmaps, all new motifs
  pdf(file = paste0(scratch, "Figures/cell_group_whole_heatmaps/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix_tmp)/44), height=2*nrow(p_cell_group_matrix_tmp)/7 )
  
  #par(pin=c(max(6,6*ncol(p_cell_group_matrix)/44), height=2*nrow(p_cell_group_matrix)/7 ))
  
  # if(newplot==FALSE){
  #   newplot=TRUE
  #   
  # }else{
  #   #plot.new()
  # }
  
  ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                             cluster_rows = TRUE,
                             cluster_columns = TRUE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             #row_names_gp = gpar(col = row_name_colors),
                             clustering_distance_rows = "pearson",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "pearson",
                             clustering_method_columns = "complete",
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
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
                               at = c(-100, 0, 100),
                               labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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
  
  result=try( draw(ht, column_title=ct_group), silent=TRUE )
  
  if( class(result)=="try-error"){
    dev.off()
    print(ct_group)
    print("error")
    pdf(file = paste0(scratch, "Figures/cell_group_whole_heatmaps/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix_tmp)/44), height=2*nrow(p_cell_group_matrix_tmp)/7 )
    ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                               cluster_rows = FALSE,
                               cluster_columns = FALSE, 
                               show_row_names = TRUE,
                               show_column_names = TRUE,
                               column_names_side = "top",
                               row_names_side = "left",
                               #row_names_gp = gpar(col = row_name_colors),
                               clustering_distance_rows = "pearson",
                               clustering_method_rows = "complete",
                               clustering_distance_columns = "pearson",
                               clustering_method_columns = "complete",
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
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
                                 at = c(-100, 0, 100),
                                 labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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
    result=try( draw(ht, column_title=ct_group), silent=TRUE )
    
  }
  
  dev.off()
  
  #Only 10 or 20 most enriched
  
  most_enriched<-list()
  most_enriched_orig<-list()
  
  #Choosing the ten most enriched for each cell type maybe not sensible, as there 
  #are many more for some cell lines, try volcano plots and predictive analysis
  
  number_of_most_enriched=10
  
  for(c in colnames(p_cell_group_matrix_new_motifs)){
    
    most_enriched[[c]]=names(sort(p_cell_group_matrix_new_motifs[,c], decreasing = TRUE)[1:number_of_most_enriched])
    most_enriched_orig[[c]]=names(sort(p_cell_group_matrix_orig_new_motifs[,c], decreasing = TRUE)[1:number_of_most_enriched])
  }
  
  
  length(unlist(most_enriched)) #70
  most_enriched=unique(unlist(most_enriched)) #37
  
  length(unlist(most_enriched_orig)) #70
  most_enriched_orig=unique(unlist(most_enriched_orig)) #44
  
  #index
  index=which(rownames(p_cell_group_matrix_orig_new_motifs) %in% unlist(most_enriched_orig))
  
  p_cell_group_matrix_new_motifs=p_cell_group_matrix_new_motifs[index,] 
  #e_cell_group_matrix=e_cell_group_matrix[index,] 
  p_cell_group_matrix_orig_new_motifs=p_cell_group_matrix_orig_new_motifs[index,] 
  
  p_cell_group_matrix_tmp=pmin(pmax(p_cell_group_matrix_new_motifs, -100), 100)
  
  #shorter_names_subset <- rownames(p_cell_group_matrix_orig) %>%
  #  map_df(~filter(shorter_names, ID == .x))
  
  
  pdf(file = paste0(scratch, "Figures/cell_group_heatmaps_new_motifs/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix_tmp)/44), height=3*nrow(p_cell_group_matrix_tmp)/7 )
  
  #par(pin=c(max(6,6*ncol(p_cell_group_matrix)/44), height=2*nrow(p_cell_group_matrix)/7 ))
  
  # if(newplot==FALSE){
  #   newplot=TRUE
  #   
  # }else{
  #   #plot.new()
  # }
  
  ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                             cluster_rows = TRUE,
                             cluster_columns = TRUE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             #row_names_gp = gpar(col = row_name_colors),
                             clustering_distance_rows = "pearson",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "pearson",
                             clustering_method_columns = "complete",
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
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
                               at = c(-100, 0, 100),
                               labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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
  
  result=try( draw(ht, column_title=ct_group), silent=TRUE )
  
  if( class(result)=="try-error"){
    dev.off()
    print(ct_group)
    print("error")
    pdf(file = paste0(scratch, "Figures/cell_group_heatmaps_new_motifs/",ct_group,".pdf"), width=max(8,8*ncol(p_cell_group_matrix_tmp)/44), height=3*nrow(p_cell_group_matrix_tmp)/7 )
    ht=ComplexHeatmap::Heatmap(p_cell_group_matrix_tmp, col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                               cluster_rows = FALSE,
                               cluster_columns = FALSE, 
                               show_row_names = TRUE,
                               show_column_names = TRUE,
                               column_names_side = "top",
                               row_names_side = "left",
                               #row_names_gp = gpar(col = row_name_colors),
                               clustering_distance_rows = "pearson",
                               clustering_method_rows = "complete",
                               clustering_distance_columns = "pearson",
                               clustering_method_columns = "complete",
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
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
                                 at = c(-100, 0, 100),
                                 labels = c("< -100 (Depleted)","0",">100 (Enriched)"),
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
    result=try( draw(ht, column_title=ct_group), silent=TRUE )
    
  }
  
  dev.off()
  
  
  
  #Interactive heatmap, does not really work here
  
  
 
  # p_cell_group_matrix_saturated=p_cell_group_matrix
  # saturated_threshold=100
  # p_cell_group_matrix_saturated[p_cell_group_matrix_saturated > saturated_threshold] <- saturated_threshold
  # p_cell_group_matrix_saturated[p_cell_group_matrix_saturated < -saturated_threshold] <- -saturated_threshold
  # 
  # 
  # 
  # metadata=representatives[match(rownames(p_cell_group_matrix_orig), representatives$ID),]
  #  
  # test=matrix(rep(metadata$ID,each=7), ncol=nrow(p_cell_group_matrix), byrow=FALSE)
  # 
  # mat=t(p_cell_group_matrix_orig) #7 x 44 (63?)
  # mat[]<-paste("Motif: ", matrix(rep(metadata$ID,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE) ,
  #              "\nFamily: ", matrix(rep(metadata$Lambert2018_families,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
  #              "\nStudy: ", matrix(rep(metadata$study,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
  #              "\n-log10 pvalue:", mat, 
  #              "\nlog2FC: ", t(e_cell_group_matrix),
  #              "\nIC: ", matrix(rep(metadata$IC,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
  #              "\nlength: ", matrix(rep(metadata$length,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
  #              "\nmatch numbers: ", matrix(rep(metadata$match_numbers,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE),
  #              "\nMOODS threshold: ", matrix(rep(metadata$threshold,7), ncol = nrow(p_cell_group_matrix), byrow=TRUE)
  # )
  
  
  # heatmaply(p_cell_group_matrix_saturated,
  #           col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  #           legend=TRUE,
  #           showticklabels = T,
  #           RowSideColors = data.frame(Family=metadata$Lambert2018_families),
  #           custom_hovertext=t(mat),
  #           row_text_angle = 0,
  #           plot_method="plotly",
  #           column_text_angle = 45,
  #           subplot_margin = 0,
  #           cellnote = NULL,
  #           cellnote_color = "auto",
  #           cellnote_textposition = "middle right",
  #           cellnote_size = 12,
  #           margins = c(NA, NA, NA, NA),
  #           scale = NULL, # Scale the data
  #           symm = FALSE,  # Make the heatmap symmetric
  #           dendrogram = "none",  # Show dendrograms on both rows and columns
  #           k_row = 2,  # Show 2 clusters in row dendrogram
  #           k_col = 2,  # Show 2 clusters in column dendrogram
  #           cluster_row = FALSE,  # Cluster rows
  #           cluster_col = FALSE,  # Cluster columns
  #           xlab = "",  # X-axis label
  #           ylab = "",  # Y-axis label
  #           main = ct_group,  # Title of the plot
  #           colorbar = TRUE,  # Show color bar
  #           key.title = "-log10 p-value \n Negative (blue) indicates depletion",
  #           colorbar_xanchor = "left",
  #           colorbar_yanchor = "middle",
  #           colorbar_xpos = 1.025,
  #           colorbar_ypos = 0.25,
  #           #colorbar_len = 0.4,
  #           colorbar_thickness = 30,
  #           fixed=FALSE,
  #           fontsize_row = 14,
  #           fontsize_col = 14,
  #           side_color_colorbar_len = 0.5,
  #           prefix="",
  #           #width = NULL,
  #           #height = NULL,
  #           file=paste0("/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys_whole/",ct_group,".html")
  # )
  
  
  
  
  
}



#Volcano plots, separate figure for each cell type group

# Sample data

options(ggrepel.max.overlaps = Inf)

# Thresholds
lfc_threshold <- 1
pvalue_threshold <- 0.01

for( ct_group in names(cell_type_groups) ){

  print(ct_group)
  #ct_group="Epithelial 1"

  #grid size
  grid_size=ceiling(sqrt(length(cell_type_groups[[ct_group]])))

  plots<-list()
  plots_new<-list()

  for(ct in cell_type_groups[[ct_group]]){
    #ct=cell_type_groups[[ct_group]][1] 
    data <- data.frame(
      name = rownames(e_matrix),
      shorter_name=shorter_names$symbol,
      log2FoldChange = e_matrix[,ct], # random log2 fold changes
      minusLog10Pvalue  = -p_matrix[,ct],
      new_motif=shorter_names$new_motif # 
    )
    
    data$pvalue <- 10^{-data$minusLog10Pvalue}
    
    
    data$significant=FALSE
    data$significant[(data$log2FoldChange > lfc_threshold | data$log2FoldChange < -lfc_threshold) & data$pvalue < pvalue_threshold]=TRUE
    
    data$new_significant=FALSE
    data$new_significant[ (data$log2FoldChange > 0.75) & data$pvalue < pvalue_threshold & data$new_motif==TRUE]=TRUE
    
    significant_motifs<- data[(data$log2FoldChange > lfc_threshold | data$log2FoldChange < -lfc_threshold) & data$pvalue < pvalue_threshold, ]
    significant_motifs_new<- data[(data$log2FoldChange > 0.75 ) & data$pvalue < pvalue_threshold & data$new_motif==TRUE, ]
    
    #non_significant_motifs<- data[-((data$log2FoldChange > lfc_threshold | data$log2FoldChange < -lfc_threshold) & data$pvalue < pvalue_threshold), ]
    
    #highlight_color <- ifelse(significant_motifs$new_motif, "red", "black")
    
    
    # plots[[ct]]=ggplot(data, aes(x=log2FoldChange, y=minusLog10Pvalue)) +
    #   geom_point(aes(color = (abs(log2FoldChange) > lfc_threshold & pvalue < pvalue_threshold)), 
    #              alpha = 0.6) +
    #   geom_point(data = subset(data, data$new_motif), color = "red", alpha = 0.6) +
    #   scale_color_manual(values = c("FALSE"="grey", "TRUE"="black")) + 
    #   geom_hline(yintercept = -log10(pvalue_threshold), linetype="dashed") +
    #   geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype="dashed") +
    #   geom_text_repel(data=subset(significant_motifs, new_motif), aes(label=shorter_name ) ,color="red")+
    #   geom_text_repel(data=subset(significant_motifs, !new_motif), aes(label=shorter_name ) ,color="black")+
    #   #geom_text(data=significant_motifs, aes(label=shorter_name), vjust=-0.5, hjust=0.5, size=3.5, check_overlap = TRUE) +
    #   theme_minimal() +
    #   theme_minimal() +
    #   labs(title=ct, x="Log2 Fold Change", y="-Log10 p-value", color="Significant") 
    
    plots[[ct]]=ggplot(data, aes(x = log2FoldChange, y = minusLog10Pvalue)) +
      
      # Plotting all points
      geom_point(aes( color = new_motif, alpha = significant)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey"))+
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype="dashed") +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype="dashed") +
      geom_text_repel(data=subset(significant_motifs, new_motif), aes(label=shorter_name ) ,color="red",box.padding = unit(0.5, "lines"))+
      geom_text_repel(data=subset(significant_motifs, !new_motif), aes(label=shorter_name ) ,color="black",box.padding = unit(0.5, "lines"))+
      theme_minimal() +
      labs(title=ct, x="Log2 Fold Change", y="-Log10 p-value",color="New motif", alpha="Significant") 
    
    
    plots_new[[ct]]=ggplot(data, aes(x = log2FoldChange, y = minusLog10Pvalue)) +
      
      # Plotting all points
      geom_point(aes( color = new_motif, alpha = new_significant)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey"))+
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
      geom_hline(yintercept = -log10(pvalue_threshold), linetype="dashed") +
      geom_vline(xintercept = c(-0.75, 0.75), linetype="dashed", col="red") +
      geom_text_repel(data=subset(significant_motifs_new, new_motif), aes(label=shorter_name ) ,color="red",box.padding = unit(0.5, "lines"))+
      #geom_text_repel(data=subset(significant_motifs, !new_motif), aes(label=shorter_name ) ,color="black",box.padding = unit(0.5, "lines"))+
      theme_minimal() +
      labs(title=ct, x="Log2 Fold Change", y="-Log10 p-value", color="New motif", alpha="Significant") 
    
    
    
    
    
  }
  
  pdf(file = paste0(scratch, "Figures/cell_group_volcanoes/",ct_group,".pdf"), width=8*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top = ct_group) )   
  dev.off()
  
  pdf(file = paste0(scratch, "Figures/cell_group_volcanoes_new_motifs/",ct_group,".pdf"), width=8*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_new, ncol = grid_size, top=ct_group))  
  dev.off()
    
}





#cluster the heatmap

#dim(p_matrix_depleted)
ed <- dist(t(p_matrix_depleted)) #Euclidean distance is the default, N times D matrix as input

# do the hierarchical clusterings
#method: the agglomeration method to be used. 
#This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#"single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

#This clusters now the cell lines based on TF binding
ehc <- hclust(ed,method="complete")
#ehc2 <- hclust(ed,method="single")
#ehc3 <- hclust(ed,method="ward.D")
#ehc4 <- hclust(ed,method="ward.D2")
#ehc5 <- hclust(ed,method="average")
#ehc6 <- hclust(ed,method="mcquitty")
#ehc7 <- hclust(ed,method="median")
#ehc8 <- hclust(ed,method="centroid")
# plot the resulting dendrogram: complete linkage method produces clear groups (2 or 3 clear clusters), 
#while single linkage starts from a single patient and goes to the most dissimilar which is not very informative here.
plot(ehc,labels=FALSE)
#plot(ehc2,labels=FALSE)
#plot(ehc3,labels=FALSE)#this looks nice, ward.D
#plot(ehc4,labels=FALSE)
#plot(ehc5,labels=FALSE)
#plot(ehc6,labels=FALSE)
#plot(ehc7,labels=FALSE)
#plot(ehc8,labels=FALSE)

## Silhouette measure plots
#Display the silhouette measure for the Clustering result based on Euclidean distance 
#for $K=2,\cdots,5$ clusters on the same plot. What is the optimal number of cluster based on the silhouette plots? 
#Hint: Search online on the meaning of silhouette plots.

#Silhouette refers to a method of interpretation and validation of consistency within clusters of data. 
#The technique provides a succinct graphical representation of how well each object has been classified.[1] 

#The silhouette value is a measure of how similar an object is to its own cluster (cohesion) 
#compared to other clusters (separation). The silhouette ranges from −1 to +1, 
#where a high value indicates that the object is well matched to its own cluster 
#and poorly matched to neighboring clusters. If most objects have a high value, 
#then the clustering configuration is appropriate. If many points have a 
#low or negative value, then the clustering configuration may have too many or too few clusters.



#We will do three things in one nested command:
#First (inner-most brackets), find the clusters using cutree(). 
#Second, calculate the silhouette measure with silhouette(). Third (outer-most brackets), make the plot.


# to demonstrate the principle, we can run the last command (for 5 clusters) in three separate commands (unhashtag the lines to see how it works)
ct <- cutree(ehc,5)
sh <- silhouette(ct,ed)
plot(sh)

plot(silhouette(cutree(ehc,2),ed))
plot(silhouette(cutree(ehc,3),ed))
plot(silhouette(cutree(ehc,4),ed))
plot(silhouette(cutree(ehc,5),ed))
plot(silhouette(cutree(ehc,6),ed))
plot(silhouette(cutree(ehc,7),ed))
plot(silhouette(cutree(ehc,8),ed))
plot(silhouette(cutree(ehc,9),ed))
plot(silhouette(cutree(ehc,10),ed))


#Display the results by a heatmap. Use the R function pheatmap for your subset of 1,500 genes with the most variance from Task 2.1. 
#Add annotations of the first four mutations found in the clinical data. 
#Use "complete" as the clustering method. Include the optimal number of clusters found in Task 2.3. 
#What is visualized in addition to the dendrogram on your heatmap? 
#  ```{r heatmap, fig.cap = 'Pretty Heatmap with Sample Annotation', warning=FALSE}
# generate an annotation data.frame to visualize the first four mutations
#anno <- data.frame(clinic[,27:30]) #176 x 90

# plot heatmap with sample annotation. Heat map produces both clustering info and gene expression info. 
#The colors correspond to the gene expression levels. Additional annotation allows to see certain genes' mutation status in patients.

min(p_matrix_depleted) #-264.5894
max(p_matrix_depleted) #1578.851


paletteLength=100
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdBu")))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))


rownames(p_matrix_depleted)=shorter_names$symbol

#CAP-selex motifs for which both TFs from the same protein family

CAP_index=which(shorter_names$experiment=="CAP-SELEX")

tmp=shorter_names[CAP_index, "Lambert2018.families"]

res = strsplit(tmp, "_")
str(res)  

res=do.call(rbind, res)

CAP_index=CAP_index[which( res[,1]==res[,2])] #120

shorter_names <- shorter_names %>% mutate(family_alt=Lambert2018.families)

shorter_names[CAP_index, "family_alt"]=res[ res[,1]==res[,2],1]

tmp=shorter_names %>% count(family_alt, sort=TRUE)
rm=unique( c( grep("_", tmp[[1]]), grep(";", tmp[[1]])) )

tmp=tmp[-rm,]
#choose only first 10
tmp=tmp[1:10,"family_alt"]
anno=shorter_names
  
#remove those family_alts that do not match tmp
length(which( !(anno$family_alt %in% tmp))) #676

anno$family_alt[!(anno$family_alt %in% tmp)]=NA

annos=data.frame(family_alt=anno$family_alt)
rownames(annos)=shorter_names$ID

# Specify colors
library("RColorBrewer")
anno_cols=c("blue", "#FF7F00", "#66A61E",
            "#E41A1C", "purple", "brown", 
            "#E7298A", "gray", "#FFE528", "cyan", "white")

names(anno_cols)=c(tmp, NA)

ann_colors = list(family_alt=anno_cols)


pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-correlation-complete.pdf"), width=150, height=15 )
pdf(file = paste0(scratch,"/Figures/top300000/rest-adult-cCREs-heatmap-correlation-complete.pdf"), width=35, height=7 )
pheatmap(p_matrix_depleted, color=color,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE,
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()

pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-correlation-complete.pdf"), width=150, height=15 )
pheatmap(t(p_matrix_depleted), color=color,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE,
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()

p_matrix_dimers=p_matrix_depleted[which(rownames(p_matrix) %in% rep_dimers),]
p_matrix_Yimeng=p_matrix_depleted[which(rownames(p_matrix) %in% rep_Yimeng),]

pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-dimers-correlation-complete.pdf"), width=100, height=20 )
pheatmap(t(p_matrix_dimers), color=color,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE, col_names_side = "top", col_hclust_side = "top",
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()

pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-Yimeng-correlation-complete.pdf"), width=50, height=20 )
pheatmap(t(p_matrix_Yimeng), color=color,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE,
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()

pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-dimers-euclidean-complete.pdf"), width=100, height=20 )
pheatmap(t(p_matrix_dimers), color=color,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE, col_names_side = "top", col_hclust_side = "top",
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()

pdf(file = paste0(scratch,"/Figures/top300000/adult-cCREs-heatmap-Yimeng-euclidean-complete.pdf"), width=50, height=20 )
pheatmap(t(p_matrix_Yimeng), color=color,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
         drop_levels=TRUE,show_rownames=TRUE,
         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
dev.off()


#pdf(file = paste0(scratch,"/Figures/adult-cCREs-heatmap-dimers-correlation-ward.pdf"), width=100, height=20 )
#pheatmap(t(p_matrix_dimers), color=color,clustering_distance_rows = "correlation",
#         clustering_distance_cols = "correlation",clustering_method="ward.D", #annotation=anno one could add protein family names, dimers vs monomers
#         drop_levels=TRUE,show_rownames=TRUE, col_names_side = "top", col_hclust_side = "top",
#         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
#dev.off()

#pdf(file = paste0(scratch,"/Figures/adult-cCREs-heatmap-Yimeng-correlation-ward.pdf"), width=50, height=20 )
#pheatmap(t(p_matrix_Yimeng), color=color,clustering_distance_rows = "correlation",
#         clustering_distance_cols = "correlation",clustering_method="ward.D", #annotation=anno one could add protein family names, dimers vs monomers
#         drop_levels=TRUE,show_rownames=TRUE,
#         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )
#dev.off()



#pdf(file = "../Figures/adult-cCREs-heatmap_euclidean_ward.pdf", width=50, height=150 )
#pheatmap(p_matrix_depleted, color=color,clustering_distance_rows = "euclidean",
#         clustering_distance_cols = "euclidean",clustering_method="ward.D", #annotation=anno one could add protein family names, dimers vs monomers
#         drop_levels=TRUE,show_rownames=TRUE, annotation_row=annos, 
#         annotation_colors=ann_colors,
#         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )

#dev.off()

#pdf(file = "../Figures/adult-cCREs-heatmap_euclidean_complete.pdf", width=50, height=150 )
#pheatmap(p_matrix_depleted, color=color,clustering_distance_rows = "euclidean",
#         clustering_distance_cols = "euclidean",clustering_method="complete", #annotation=anno one could add protein family names, dimers vs monomers
#         drop_levels=TRUE,show_rownames=TRUE, annotation_row=annos, 
#         annotation_colors=ann_colors,
#         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )

#dev.off()

#pdf(file = "../Figures/adult-cCREs-heatmap_correlation_ward.pdf", width=50, height=150 )
#pheatmap(p_matrix_depleted, color=color,clustering_distance_rows = "correlation",
#         clustering_distance_cols = "euclidean",clustering_method="ward.D", #annotation=anno one could add protein family names, dimers vs monomers
#         drop_levels=TRUE,show_rownames=TRUE,annotation_row=annos, 
#         annotation_colors=ann_colors,
#         show_colnames=TRUE,legend=TRUE, breaks= myBreaks) #, cutree_rows=6, )

#dev.off()



#cutree_cols=2)













library("RColorBrewer")
paletteLength=100
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdBu")))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))


#Experimental
library(d3heatmap)
#install.packages("d3heatmap")
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("talgalili/d3heatmap")

# Generate sample data
data <- matrix(rnorm(100), nrow = 10)

# Create interactive heatmap
d3heatmap(t(epithelial), 
          scale = "none",
          dendrogram = "none",
          #colors = ength=100,
          #color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(paletteLength),
         # k_row = 2,
          rng=c(100,-100),
          k_col = 2,
          cluster_row = FALSE,
          cluster_col = FALSE,
          key = TRUE,
          main = "Interactive Heatmap with d3heatmap",
         color = "RdBu",  # Specify color palette
          na_color = "white",  # Specify color for missing values
          show_dendro = TRUE,  # Show dendrograms
          hoverCol = "yellow",  # Specify hover color
          hoverRow = "yellow",  # Specify hover color
          dendro_height = 300,  # Specify dendrogram height
          dendro_width = 100,  # Specify dendrogram width
         breaks=myBreaks
)




d3heatmap(
  t(epithelial),
  main = "Epithelial cells",
  col=colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
  #col=hcl.colors(12, "RdBu", rev = TRUE),
  width = NULL,
  height = NULL,
  show_grid = FALSE,
  anim_duration = 500,
  rng = c(-100,100),
  symm = FALSE,
  Rowv = TRUE,
  distfun = NULL,
  hclustfun = NULL,
  dendrogram = "none", #c("both", "row", "column", "none"),
  reorderfun = function(d, w) reorder(d, w),
  k_row = NULL,
  k_col = NULL,
  revC = FALSE,
  scale = "none",# c("none", "row", "column"),
  scale.by.range = FALSE,
  na.rm = FALSE,
  na.color = "#777777",
  na.value = NA,
  keysize = 1,
  key.title = "- log10 p-value",
  key.location = "fl", #c("fl", "br", "tr", "tl", "bl"), float(follows the location of axis labels key = TRUE,) bottom_right top_right top_left bottom_left
  density.info = "histogram",c("histogram", "none"),
  denscol = NULL,
  digits = 3L,
  cellnote = NULL,
  cellnote_scale = FALSE,
  cellnote_row = NULL,
  cellnote_col = NULL,
  cellnote_val = "Value",
  brush_color = "#0000FF",
  print.values = FALSE,
  notecol = "#222222",
  notecex = 1,
  #cex.note,
  theme = NULL,

  #col = "RdBu",
  breaks = c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
             seq(100/paletteLength, max(100), length.out=floor(paletteLength/2))),
  symbreaks = TRUE,
  #labRow = rownames(x),
  #labCol = colnames(x),
  labColSize = 80,
  labRowSize = 120,
  cexCol = NULL,
  cexRow = NULL,
  srtCol = 60,
  sideCol = 4,
  sideRow = 2,
  #xaxis.location,
  #yaxis.location,
  xlab = "Motif",
  ylab = "Cell type",
  xaxis_title_font_size = 14,
  yaxis_title_font_size = 14,
  ColSideColors = NULL,
  RowSideColors = NULL,
  RowColorsPalette = NULL,
  ColColorsPalette = NULL,
  RowColorsNames = NULL,
  ColColorsNames = NULL,
  key = TRUE,
)




# Generate sample data
#data <- matrix(rnorm(100), nrow = 10)

gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0, limits = c(-100, 100)
)





# Generate sample data
data <- matrix(rnorm(100), nrow = 10)

# Define color palette using RColorBrewer
palette <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)

# Create heatmap with custom color palette
heatmaply(data, 
          scale = NULL,
          symm = TRUE,
          dendrogram = "both",
          k_row = 2,
          k_col = 2,
          cluster_row = TRUE,
          cluster_col = TRUE,
          xlab = "Columns",
          ylab = "Rows",
          main = "Interactive Heatmap with Custom Colors",
          colorbar = TRUE,
          breaks=seq(-2, 2, by = 0.5), 
          col = palette  # Specify custom color palette
)

       
# set up panel of 2x2 plots:
par(mfrow=c(2,2))

# plots for 2, 3, and 4 clusters:



#Which are representative motifs

representatives=read.table("../../../motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]
rep_motifs2=rep_motifs

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

length(unique(rep_motifs))
#1062



cCREs_list<- readRDS(file = paste0( "../RData/cCREs_list.Rds") ) #
length(cCREs_list)

all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs)


N=length(all_cCREs) #all possible cCREs #435142



library(ggplot2)
# Basic histogram

df=data.frame(cCREs_nro= unlist(lapply(cCREs_list, length)))

max(df$cCREs_nro) #48391
min(df$cCREs_nro) #846

#ggplot(df, aes(x=_nro)) + geom_histogram()
# Change the width of bins
pdf(file = "../Figures/cCREs-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df, aes(x=cCREs_nro)) + geom_histogram(binwidth=1000)
dev.off()

#Width of cCREs, these are all the same 400 bp?

widths=lapply(cCREs_list, width)
lapply(widths, table)

match_numbers=lapply(representative_motif_matches_Vierstra, length)

df2=data.frame(match_nro=unlist(match_numbers))

max(df2$match_nro) #4063189
min(df2$match_nro) #21108


# Change the width of bins
pdf(file = "../Figures/Vierstra-representative-motif-hit-nros.pdf" )
#width, height, onefile, family, title, fonts, version,
#paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#useDingbats, useKerning, fillOddEven, compress)
ggplot(df2, aes(x=match_nro)) + geom_histogram(binwidth=10000)
dev.off()

mean(df2$match_nro) #307578
median(df2$match_nro) #247368

#Are the low hitters longer motifs/dimers and super hitters short

#select only 1062 representatives
representatives_match_numbers=representatives[which(representatives$new_representative=="YES"),]
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives_match_numbers$filename,"/"))[,5])

table(test==rownames(df2) ) #the motif names are in the same order

representatives_match_numbers$match_numbers=df2$match_nro

monomers=c("monomeric", "monomer or dimer", "monomer", "monomeric or dimeric")

multimers=c("dimeric", "putative multimer", "putatively multimeric", "dimer of dimers", "trimeric")


#both monomeric

representatives_match_numbers$motif_type=rep("monomer", nrow(representatives_match_numbers))

cap_selex=representatives_match_numbers$experiment=="CAP-SELEX"

mm=representatives_match_numbers$type %in% multimers

representatives_match_numbers[mm | cap_selex, "motif_type"]="multimer"

# Overlaid histograms
ggplot(representatives_match_numbers, aes(x=match_numbers, color=motif_type)) + geom_histogram(fill="white",binwidth=10000)+ xlim(0, 1e6)

#There seems to be not much difference between the distributions between monomers and multimers

#Means and medians of monomers and multimers

representatives_match_numbers %>%
  group_by(motif_type) %>%
  dplyr::summarize(Mean = mean(match_numbers, na.rm=TRUE), Median = median(match_numbers, na.rm=TRUE))

#lengths of the motifs

representatives_match_numbers$motif_lengths=lapply(representative_motif_matches_Vierstra, function(x) unique(width(x)) )

#count as the function of motif length

ggplot(representatives_match_numbers, aes(x=as.numeric(motif_lengths), y=as.numeric(match_numbers))) + geom_point()+ geom_smooth()

+scale_x_continuous()+scale_y_continuous()


#how many have more than half million

length(which(representatives_match_numbers$match_numbers>500000)) #128


#what are these motifs
super_hitters=representatives_match_numbers[representatives_match_numbers$match_numbers> 500000,]
  
  
#What are those with low number of hits
low_hitters=representatives_match_numbers[representatives_match_numbers$match_numbers< 50000,]


#distribution of scores for high hitters and low hitters

#the highest hitter
arrange(super_hitters, desc(match_numbers))$ID[1]
#why there is motif ID, remove
hist(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score)

#How many if using score threshold 5
table(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score>=5)

hist(representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score[representative_motif_matches_Vierstra[[ arrange(super_hitters, desc(match_numbers))$ID[1] ]]$score>=5])

#the lowest hitter
arrange(low_hitters, match_numbers)$ID[1]
hist(representative_motif_matches_Vierstra[[ arrange(low_hitters, match_numbers)$ID[1] ]]$score)


for(hits in motifs){
  print(hits)
  #hits=motifs[1]

 
  #m: the size of all cCREs that have a motif match#
  fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], all_cCREs, type="within", ignore.strand=TRUE)
  
  m=length(unique(fol@to)) #This is now number of all cCREs that have at least one motif hit #9857
  
  for(cell_type in names(cCREs_list)){
    #cell_type=names(cCREs_list)[1]
    #print(cell_type)
    
    #findOverlaps(query, subject)
    fol=findOverlaps(representative_motif_matches_Vierstra[[hits]], cCREs_list[[cell_type]], type="within", ignore.strand=TRUE)
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


saveRDS(p_matrix,  file = paste0( "../RData/Vierstra_p_matrix_human_cell_type_restricted.Rds")) #
saveRDS(e_matrix,  file = paste0( "../RData/Vierstra_e_matrix_human_cell_type_restricted.Rds")) #





