library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library("tidyverse")
#library("pheatmap")
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


#.libPaths(c("/projappl/project_2006203/project_rpackages_4.2.1"))
.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new")
#library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")


setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
data_path="/scratch/project_2006203/TFBS/"

#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

#p_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs.Rds")) # #The TFs are in alphabetical order
e_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs.Rds")) #


scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

#source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")
#source("../code/heatmap_motor.R")

#Cell type names mapping

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Adult_Celltypes_mapping <- read_delim("../CATLAS/Adult_Celltypes_mapping.csv", 
#                                       delim = ";", escape_double = FALSE, trim_ws = TRUE)


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





colnames(e_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(e_matrix), Adult_Celltypes_mapping$celltype) ]



#Logistic regression coefficients

representatives$new_motif=FALSE

representatives$new_motif[representatives$study=="fromYimeng"]=TRUE 

regression_coefs=matrix(0, nrow=length(which(representatives$new_representative=="YES")),
                        ncol=sum(sapply(cell_type_groups, length))) #motifs as rows, columns as cell types
rownames(regression_coefs)=representatives$ID[which(representatives$new_representative=="YES")]
colnames(regression_coefs)=colnames(e_matrix)

library(glmnet)

for( cell_type_group in names(cell_type_groups) ){

  #cell_type_group_nro=1
  #cell_type_group=names(cell_type_groups)[cell_type_group_nro]

  print(cell_type_group)
  load(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/logistic_regression_processed_",cell_type_group,".RData"))

  for(cell_type in names(cvfit_final_list_maxscore) ){
    print(cell_type)
    #cell_type=names(cvfit_final_list_maxscore)[1]
    regression_coef=data.frame(name=rownames(coef(cvfit_final_list_maxscore[[cell_type]], s = "lambda.1se"))[-1],
                           coef=as.matrix(coef(cvfit_final_list_maxscore[[cell_type]], s = "lambda.1se"))[-1] )
    regression_coefs[,cell_type]=regression_coef$coef
  }
}

saveRDS(regression_coefs, paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/all_regression_coefs_maxscore.RDS"))


