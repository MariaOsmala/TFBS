library("rtracklayer")
library("dbplyr")
library("dplyr")
library("tidyverse")
library("ComplexHeatmap")
library("cluster")
library("readr")
library("grid")
library("RColorBrewer")
library("ggrepel") #
library("circlize")
library(gridExtra)
library(grid) # for textGrob



.libPaths(c("/projappl/project_2006203/project_rpackages_4.3.0/"))
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
#setwd("~/projects/TFBS/ATAC-seq-peaks/RProjects/")
data_path="/scratch/project_2006203/TFBS/"
#data_path="~/projects/TFBS/"


representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3933

bHLHhomeo_all=representatives[which(representatives$Lambert2018_families %in% c("bHLH_Homeodomain","Homeodomain_bHLH")),] #59
rep_motifs=bHLHhomeo_all$ID

p_matrix=readRDS(file = paste0(data_path,"ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs_version2.2_bHLH_homeo.Rds")) #
e_matrix=readRDS(file = paste0(data_path, "/ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs_version2.2_bHLH_homeo.Rds")) #




source("../code/heatmap_motor.R")

#Cell type names mapping

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#Adult_Celltypes_mapping <- read_delim("../CATLAS/Adult_Celltypes_mapping.csv", 
#                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)


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

 
#version 2.2

# The new motifs, select only representatives bHLH-homeo CAP-SELEX motifs
# that represent the new bHLH-Homeodomain spacing motifs

# bHLHhomeo_all contain 59 motifs

# bHLH_homeo_new contain 40 motifs
# 33 spacing (13 representatives)
# 7 composite (4 representatives)
# 17 representatives
# 23 non-representatives

bHLHhomeo_new=bHLHhomeo_all %>% filter(study=="fromYimeng") #40
bHLHhomeo=bHLHhomeo_new$ID 
table(bHLHhomeo_new$new_representative) #17 representatives
bHLHhomeo_new %>% filter(new_representative=="YES") %>% count(type)# 13 spacing, 4 composite

#Select only the spacing motifs
bHLHhomeo=bHLHhomeo_new %>% filter(type=="spacing")  # 33 spacing, 13 representatives, 20 non representatives

write_delim(bHLHhomeo_all %>% filter(study=="fromYimeng" & type=="spacing") %>% select(ID), file="bHLH_homeodomain_new_spacing_all.tsv", col_names = FALSE) #33
write_delim(bHLHhomeo %>% filter(new_representative=="YES") %>% select(ID), file="bHLH_homeodomain_new_spacing_representatives.tsv", col_names = FALSE) #13

#What are the old CAP-SELEX motifs that represent the new non-representatives spacing motifs?
represented_by_representatives <- readRDS("/projappl/project_2006203/TFBS/RProjects/TFBS/RData/represented_by_representatives_version2.2.RDS")
non_reps=bHLHhomeo %>% filter(new_representative=="NO") %>% pull(ID) #20

indices=list()
for(non_rep in non_reps){
  indices[[non_rep]] <- names(which(sapply(represented_by_representatives, function(x) any(grepl(non_rep, x))))  )
}

# 1. Find the maximum length of the list elements
max_length <- max(sapply(indices, length))

# 2. Pad shorter elements with NA to match the longest element
padded_list <- lapply(indices, function(x) {
  length(x) <- max_length
  x
})

# 3. Combine the list elements into a matrix/data frame
indices_df <- as.data.frame(do.call(rbind, padded_list))
indices_df$ID=rownames(indices_df)
rownames(indices_df)=NULL
names(indices_df)=c("rep1", "rep2", "ID")
indices_df=indices_df[c("ID","rep1", "rep2")]

reps=c(as.vector(indices_df$rep1), as.vector(indices_df$rep2))
reps=unique(reps[!is.na(reps)]) #7

#Are all these already in the set bHLHhomeo
reps %in% bHLHhomeo$ID

#Are they in the bHLHhomeo all
reps %in% bHLHhomeo_all$ID

#which of these representatives are in the new set
representatives %>% filter(ID %in% reps) %>% select(study)
#which are in the old set
representatives %>% filter(ID %in% reps & study!="fromYimeng") %>% select(ID)

#Remove the HT-SELEX motif from Yin 2017
reps=reps[-which(reps=="NEUROG1_HT-SELEX_TTGCCA40NGTG_KAE_RNCATATGNY_1_2")]

#Are all these already in the set bHLHhomeo
reps %in% bHLHhomeo$ID

#Are they in the bHLHhomeo all
missing_reps=reps[-which(reps %in% bHLHhomeo$ID)]

missing_reps %in% bHLHhomeo_all$ID

bHLHhomeo_reps=bHLHhomeo %>% filter(new_representative=="YES") #13

bHLHhomeo_reps=rbind(bHLHhomeo_reps, bHLHhomeo_all %>% filter(ID %in% missing_reps) ) #15

bHLHhomeo_reps$ID

#Remove
#Human TALE (Three Amino acid Loop Extension) homeobox genes for an "atypical" homeodomain
#"BHLHE22_TGIF2LX_TAACAG40NGGC_YVII_NCATATGNNNNNTGTCAN_m1_c3b0" TGIF2LX is TALE? Homeobox; TALE-type 
#"ATOH1_HOXC10_TTCGGG40NGAG_YAEIII_NYMRTAAANNNNNCATATGN_m1_c3b0" 
#"HOXD10_TFEC_TAGCTG40NCTA_YAGIIII_MRTAAANNNNNNNCACGTG_m1_c3b0"


remove_ind=which(bHLHhomeo_reps$ID %in% c("BHLHE22_TGIF2LX_TAACAG40NGGC_YVII_NCATATGNNNNNTGTCAN_m1_c3b0",
"ATOH1_HOXC10_TTCGGG40NGAG_YAEIII_NYMRTAAANNNNNCATATGN_m1_c3b0", 
"HOXD10_TFEC_TAGCTG40NCTA_YAGIIII_MRTAAANNNNNNNCACGTG_m1_c3b0"))

bHLHhomeo_reps=bHLHhomeo_reps[-remove_ind,] #12


final_matrix=p_matrix_depleted[ which(rownames(p_matrix_depleted) %in% bHLHhomeo_reps$ID), ]

#remove uninteresting cells
final_matrix=final_matrix[,-which( colSums( ( abs(final_matrix) < 50  | (final_matrix < 0) ) )== nrow(final_matrix) )]




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
row_name_info=as.data.frame(row_name_info)
row_name_info$ID=rownames(final_matrix)

row_name_info[1,3:7]=row_name_info[1,4:8]
row_name_info[1,8]=NA

# Regular expression for a sequence of 'N's in the middle

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

rownames(final_matrix)=apply(row_name_info[,c(1,2,10)],1,paste0,collapse="_")

#Rotate column labels to that they can be read from up to down
#Add protein family info

fontsize=12
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

ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                        #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                        col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
                        cluster_rows = TRUE,
                        cluster_columns = TRUE, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_side = "top",
                        row_names_side = "left",
                        column_names_rot=45, #315,
                        #row_names_gp = gpar(col = row_name_colors),
                        #row_title_gp=gpar(fontface=row_name_fontfaces ),
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


ht_list=draw(ht_list)

dev.off()

r.dend <- row_dend(ht_list)  #Extract row dendrogram
rcl.list <- row_order(ht_list)  #Extract clusters (output is a list)

rownames(final_matrix[rcl.list,])

row_order=c("BHLHA15_EN2_TAATTA3NCATATG",  #3

"BHLHB8_ISX_NTAATTR4NCAKMTGN",  
"ATOH1_HOXB5_NTAATTA4NCAKMTGN", #4
"HOXB2_TCF3_NCACCTG5NMATTA",  
"BARHL1_TCF12_NTAATTG4NCACCTGN",

"NHLH1_EVX2_NMNCAKSTG5NTAATTAN", #5
"OLIG2_EVX2_NTAATTA5NCATATGN", 

"OLIG2_HOXB5_NCATATG6NYMATTAN",  #6
"ATOH1_EVX2_NTAATTA6NCAKMTGN",  
"NHLH1_EVX2_NMNCAKSTG6NTAATTAN",
"CLOCK_ISL2_NCAMTTA6NCACGTGN", 

"NHLH1_HOXD8_NCASCTG8NYAATTRN") #8


row_split = rep("6", length(row_order))
row_split[1] = "3"
row_split[2:5]="4"
row_split[6:7]="5"
row_split[12]="8"

row_split=factor(row_split, levels=c("3", "4", "5", "6", "8"))

final_matrix=final_matrix[row_order,]

#Add the protein family and cell_type_group info

row_name_info$new_name_full=apply(row_name_info[,c(1,2,10)],1,paste0,collapse="_")

#Stromal first, then Islet, Then Brain glia, Then Brain neuron  
column_order=c("Fibroblast (Epithelial)", "Peripheral Nerve Stromal", "Fibroblast (General)" ,   "Pancreatic Beta Cell 1", "Pancreatic Beta Cell 2", 
                "Pancreatic Delta,Gamma cell","Pancreatic Acinar Cell",
               "Plasma Cell",  
               "Astrocyte 1", "Astrocyte 2",  "Glutamatergic Neuron 1", "Glutamatergic Neuron 2",  "GABAergic Neuron 1", "GABAergic Neuron 2")     


final_matrix=final_matrix[, column_order]

#CEll type groups
ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]

#One needs to order the row_name_info also
row_name_info=row_name_info[match( row_order,row_name_info$new_name_full),]

#Set symbols as rownames

rownames(final_matrix)=paste0(row_name_info$V1, "_",row_name_info$V2)


#Highlight new representative motifs by bold font
bold=representatives$study[match(row_name_info$ID, representatives$ID)]
row_name_fontfaces <- ifelse(bold %in% c("fromYimeng"), "bold", "plain")


# Generate a color palette
num_colors <- length(unique(ct_group_names))
palette_colors <- brewer.pal(num_colors, "Set2")  # Or any other palette

# Create a named vector of colors
annotation_colors <- setNames(palette_colors, unique(ct_group_names))

image_png = paste0("../../PWMs_final_version2.2/Logos_puhti/Logos_final2/png/prob/", row_name_info$ID, ".png")
image_png=paste0("/scratch/project_2006203/TFBS/PWMs_final_version2.2/Logos_final2/png/prob/", row_name_info$ID, ".png")
#image_pdf = paste0("../../PWMs_final_version2.2/Logos_mac/Logos_final2/pdf/prob/", row_name_info$ID, ".pdf")
image_pdf = paste0("../../PWMs_final_version2.2/Logos_puhti/Logos_final2/pdf/prob/", row_name_info$ID, ".pdf")
image_pdf= paste0("/scratch/project_2006203/TFBS/PWMs_final_version2.2/Logos_final2/pdf/prob/", row_name_info$ID, ".png")

#Draw also reverse complements
image_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2_Logos/svg/", row_name_info$ID, ".svg")

image_revcomp_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2_rev_comp_Logos//svg/", row_name_info$ID, ".svg")

image_svg[c(4,6,8,10,12)]=image_revcomp_svg[c(4,6,8,10,12)]


# we only draw the image annotation for PNG images, while the others are the same
#install.packages("grImport2")
#install.packages("rsvg")

ha = HeatmapAnnotation(Logo = anno_image(image_svg, height= unit(8, "cm"), 
                                         border = FALSE,
                                         width= unit(8, "cm"),
                                         space = unit(0, "mm")), which="row") #, height=unit(0.618,"cm")) #, 



pdf(file = paste0(scratch, "Figures/bHLH_homeo_representatives.pdf"), width=14, 
    height=12)

ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                                #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                                col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
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
                                height = nrow(final_matrix)*unit(0.8,"cm"), 
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
                                  at = c(-20, 0, 100),
                                  labels = c("-20", "0",">100 (Enriched)"),
                                  labels_gp = gpar(fontsize = fontsize),
                                  #labels_rot = 0,
                                  #nrow = NULL,
                                  #ncol = 1,
                                  #by_row = FALSE,
                                  legend_gp = gpar(fontsize = fontsize),
                                  #legend_height = NULL,
                                  #legend_width = NULL,
                                  legend_direction = c("horizontal")
                                  # break_dist = NULL,
                                  #graphics = NULL,
                                  #param = NULL
                                )
)

draw(ht_list+ha, row_title="spacing", 
     row_title_gp=gpar(fontface="bold"), 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.off()

