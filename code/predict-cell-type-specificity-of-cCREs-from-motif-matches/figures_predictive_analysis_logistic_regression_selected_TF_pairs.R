
library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")


#library("GRanges")
library("readr")
library("caret")
library("glmnet")
library("reshape2")
library(ggrepel) #
library("gridExtra")
library("openxlsx")

#.libPaths("/projappl/project_2006203/project_rpackages_4.2.1")
#.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")
#install.packages("ggbreak")
library(ggbreak)

setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
data_path="/scratch/project_2006203/TFBS/"

scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

representatives$new_motif=FALSE

representatives$new_motif[representatives$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420")]=TRUE   

length(unique(rep_motifs))

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"

cell_type_groups<-list()
cell_type_groups[["Epithelial 1"]]=c("Parietal Cell",
                                     "Chief Cell",
                                     "Pancreatic Acinar Cell",
                                     "Ductal Cell (Pancreatic)",
                                     "Gastric Neuroendocrine Cell",
                                     "Foveolar Cell",
                                     "Hepatocyte")
cell_type_groups[["Stromal"]]=c("Mesothelial Cell",
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
                                "Smooth Muscle (Esophageal Muscularis) 3",
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


regression_coefs=readRDS( paste0( data_path, "ATAC-seq-peaks/RData/all_regression_coefs_maxscore.RDS"))
scaled_regression_coefs=scale(regression_coefs,center=TRUE, scale=TRUE)


#Cell type groups of interest
ct_groups=c("Stromal", "Immune Lymp", "Immune Myelo", "Endothelial",  "Brain glia")               
#cell_type_group_nro <- as.numeric(commandArgs(trailingOnly = TRUE)) #1:15


#Motifs of interest
family_pair=c("Ets_Runt", "Runt_Ets","Forkhead_Ets", "Ets_ForkHead") 
family_pair_metadata=representatives[which(representatives$Lambert2018_families %in% family_pair & representatives$new_representative=="YES"),]
family_pair_metadata=family_pair_metadata[-which(family_pair_metadata$ID=="FOXO1_ETV7_CAP-SELEX_TTGGTG40NTAC_AS_RWMAACAGGNNNNNNTTCCNN_1_2"),] #Remove strange composite

row_name_info=data.frame(ID=family_pair_metadata$ID)

#Shorter names for motifs that contain the seed

#Convert the motif names (row labels) name + 5th column, Ns shortened
my_list=strsplit(family_pair_metadata %>% filter(type=="pfm_spacing_new" ) %>% select(ID) %>%pull(ID), "_")
# Find the length of the longest vector
max_length <- max(sapply(my_list, length))

# Pad shorter vectors with NA
padded_list <- lapply(my_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert the list to a data frame
row_name_info[which(family_pair_metadata$type=="pfm_spacing_new") , 2:6] <- do.call(rbind, padded_list)[,-seq(6,10,1)]



#Do this only for the spacing motifs
# Regular expression for a sequence of 'N's in the middle
pattern <- "(?<!^)(NN+)(?!$)"

# Find positions of the pattern
matches <- gregexpr(pattern, row_name_info[which(family_pair_metadata$type=="pfm_spacing_new"),6], perl = TRUE)
matches_info=as.data.frame(cbind(row_name_info[which(family_pair_metadata$type=="pfm_spacing_new"),5],do.call(rbind, lapply(matches,attributes))))

replacement=paste0(unlist(matches_info$capture.length) ,"N")

matches_info$new_name=""

for(i in 1:nrow(matches_info)){
  #print(i)
  matches_info$new_name[i]=gsub(pattern, replacement[i], row_name_info[which(family_pair_metadata$type=="pfm_spacing_new")[i],6], perl = TRUE)
}

row_name_info$new_name=""
row_name_info$new_name[which(family_pair_metadata$type=="pfm_spacing_new") ]=matches_info$new_name

#COmposite motifs

#Convert the motif names (row labels) name + 5th column, Ns shortened


row_name_info$new_name[which(family_pair_metadata$type %in% c("pfm_composite_new", "")) ]=family_pair_metadata %>% filter(type %in% c("pfm_composite_new", "") ) %>% select(seed) %>%pull(seed)
row_name_info$symbol=family_pair_metadata$symbol
row_name_info$final_new_name=paste0(row_name_info$symbol,"_", row_name_info$new_name)

family_pair_metadata$ID
row_name_info$ID

row_ind=match(row_name_info$ID, rownames(regression_coefs))

rownames(regression_coefs)[row_ind]=row_name_info$final_new_name




for(ct_group in ct_groups){
print(ct_group)
#ct_group=ct_groups[1]
  plots=list()

for(cell_type in cell_type_groups[[ct_group]] ){
    #cell_type=cell_type_groups[[ct_group]][1]   
    
    regression_coef=data.frame(coef=regression_coefs[,cell_type], row.names = NULL)
    regression_coef$name=rownames(regression_coefs)
    # 
    regression_coef=regression_coef[order(regression_coef$coef, decreasing =TRUE),]
    # 
    regression_coef$order=1:nrow(regression_coef)
    
    regression_coef$highlight=regression_coef$name %in% row_name_info$final_new_name
    
  
        
    # 
    plots[[cell_type]]= ggplot(regression_coef)+geom_point(data = subset(regression_coef, !highlight), aes(x = order, y = coef), color = "darkgrey") +
      # Plot red points on top
      geom_point(data = subset(regression_coef, highlight), aes(x = order, y = coef), color = "red") +
      
      geom_text_repel(data=subset(regression_coef,highlight), aes(x = order, y = coef,label=name) ,color="red",
                      #box.padding = unit(0.5, "lines"),point.padding = unit(0.5, "lines"), 
                      #force=2,
                      box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.5, "lines"),
                      #segment.padding = unit(0.3, "lines"),
                      #nudge_x = 0.05,
                      #nudge_y = 0.05,
                      #size = 3
                      #max.overlaps = 5
                      
                      )+
      
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient")+
      theme_minimal()
    }
  
  #grid size
  grid_size=ceiling(sqrt(length(cell_type_groups[[ct_group]])))
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs_selected_motifs/","Ets_Runt_Forkhead_",ct_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top=ct_group))  
  dev.off()
  
  }

