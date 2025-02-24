#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#.libPaths("/projappl/project_2006203/project_rpackages_4.2.1")


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

#source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")
source("../code/heatmap_motor.R")

#short names for TFs

#Cell type names mapping

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

Adult_Celltypes_mapping <- read_delim("../CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)


#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"



scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

data_path="/scratch/project_2006203/TFBS/"


#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

p_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs.Rds")) # the motif names are in alphabetical order
e_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs.Rds")) #

#sort the representatives table according to the rownames of p_matrix
representatives <- representatives[order(match(representatives$ID, rownames(p_matrix))), ]


test=representatives$ID
rep_motifs=test[which(representatives$new_representative=="YES")]
length(which(rep_motifs %in% rownames(p_matrix))) #1031

length(unique(rep_motifs)) #1031

shorter_names=representatives %>% filter(test %in% rownames(p_matrix)) %>% select(symbol, ID, experiment, 
                                              new_representative, Lambert2018_families, type)
#table(shorter_names)

shorter_names$new_motif=FALSE
shorter_names$new_motif[which(shorter_names$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420") )]=TRUE #245


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
colnames(p_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(p_matrix), Adult_Celltypes_mapping$celltype) ]



#which are enriched 

infinites=which(is.infinite(e_matrix)) #226
e_matrix[infinites]=NA
enriched=which(e_matrix>0) #1062*111=117882


enrichments=e_matrix[enriched] #-Inf to 2.7
min(enrichments, na.rm=TRUE) #2.184901e-06
max(enrichments, na.rm=TRUE) #4.473739

#which are depleted 
depleted=which(e_matrix<0)
depletions=e_matrix[depleted]
min(depletions, na.rm=TRUE) #-3.0469
max(depletions, na.rm=TRUE) #-1.704847e-06

#Convert the depleted -log_10 p-values to negative
str(depleted)

p_matrix_depleted=-p_matrix #1031 x 11
p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]


#Volcano plots, separate figure for each cell type group

# Sample data
library("gridExtra")
options(ggrepel.max.overlaps = Inf)

# Thresholds
lfc_threshold <- 1
pvalue_threshold <- 0.01


family_pair=c("Ets_Runt", "Runt_Ets","Forkhead_Ets", "Ets_ForkHead") 
family_pair_metadata=representatives[which(representatives$Lambert2018_families %in% family_pair & representatives$new_representative=="YES"),]

family_pair_metadata=family_pair_metadata[-which(family_pair_metadata$ID=="FOXO1_ETV7_CAP-SELEX_TTGGTG40NTAC_AS_RWMAACAGGNNNNNNTTCCNN_1_2"),] #Remove strange composite





family_pair_metadata$type

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
#pattern <- "(?<=.)(N+)(?=.+$)"
#pattern <- "(?<!^)(N+)(?!$)"
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

row_ind=match(row_name_info$ID, rownames(e_matrix))
rownames(e_matrix)[row_ind]=row_name_info$final_new_name

representatives$new_motif=FALSE

representatives$new_motif[representatives$type %in% c("pfm_composite_new", "pfm_spacing_new", "20230420")]=TRUE 

ct_groups=c("Stromal", "Immune Lymp", "Immune Myelo", "Endothelial",  "Brain glia")               


for( ct_group in ct_groups ){

  print(ct_group)
  #ct_group="Endothelial"

  #grid size
  grid_size=ceiling(sqrt(length(cell_type_groups[[ct_group]])))

  plots<-list()


  for(ct in cell_type_groups[[ct_group]]){
    #ct=cell_type_groups[[ct_group]][1] 
    data <- data.frame(
      name = rownames(e_matrix),
      shorter_name=shorter_names$symbol,
      log2FoldChange = e_matrix[,ct], # random log2 fold changes
      minusLog10Pvalue  = -p_matrix[,ct],
      new_motif=representatives$new_motif # 
    )
    
    data$pvalue <- 10^{-data$minusLog10Pvalue}
    
    data$highlight=data$name %in% row_name_info$final_new_name
    
        
    plots[[ct]]=ggplot(data) +
      # Plot grey points first
      geom_point(data = subset(data, !highlight), aes(x = log2FoldChange, y = minusLog10Pvalue), color = "darkgrey") +
      # Plot red points on top
      geom_point(data = subset(data, highlight), aes(x = log2FoldChange, y = minusLog10Pvalue), color = "red") +
      #geom_hline(yintercept = -log10(pvalue_threshold), linetype="dashed") +
      #geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype="dashed") +
      geom_text_repel(data=subset(data,highlight), aes(x = log2FoldChange, y = minusLog10Pvalue,label=name) ,color="red",box.padding = unit(0.5, "lines"))+
      theme_minimal() +
      labs(title=ct, x="Log2 Fold Change", y="-Log10 p-value") 
    
    
    
    
    
    
  }
  
  pdf(file = paste0(scratch, "Figures/cell_group_volcanoes/Ets_Runt_Forkhead_",ct_group,".pdf"), width=8*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top = ct_group) )   
  dev.off()
  
    
}






