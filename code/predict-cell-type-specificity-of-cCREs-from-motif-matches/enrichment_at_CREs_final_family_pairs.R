library("rtracklayer")
library("dbplyr")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("pheatmap")
library("ComplexHeatmap")
library("cluster")
library("readr")
library("grid")
library("RColorBrewer")

library("ggrepel") #
library("circlize")
library(gridExtra)
library(grid) # for textGrob

.libPaths("/projappl/project_2006472/project_rpackages_4.3.0_new")
#library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")


setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
data_path="/scratch/project_2006203/TFBS/"

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3294

#p_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_p_matrix_human_cell_type_restricted_all_motifs.Rds")) # #The TFs are in alphabetical order
#e_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/final_e_matrix_human_cell_type_restricted_all_motifs.Rds")) #

scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

#source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")
source("../code/heatmap_motor.R")

#Cell type names mapping

source("../code/cell_type_names.R")

length(cell_type_groups) #15, in the paper they mention 17


#colnames(e_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(e_matrix), Adult_Celltypes_mapping$celltype) ]
#colnames(p_matrix)=Adult_Celltypes_mapping$`Cell type`[ match(colnames(p_matrix), Adult_Celltypes_mapping$celltype) ]

#which are enriched 

# infinites=which(is.infinite(e_matrix)) #226
# e_matrix[infinites]=NA
# 
# enriched=which(e_matrix>0) 
# 
# enrichments=e_matrix[enriched] #-Inf to 2.7
# min(enrichments, na.rm=TRUE) #2.184901e-06
# max(enrichments, na.rm=TRUE) #4.473739

#which are depleted 
# depleted=which(e_matrix<0)
# depletions=e_matrix[depleted]
# min(depletions, na.rm=TRUE) #-3.0469
# max(depletions, na.rm=TRUE) #-1.704847e-06
# 
# p_matrix_depleted=-p_matrix #1031 x 11
# p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]


#family_pair=c("Ets_Runt", "Runt_Ets") #These are composite motifs 3

family_pair=c("Forkhead_Ets", "Ets_ForkHead") #These are composite motifs 4

family_pair_metadata=representatives[which(representatives$Lambert2018_families %in% family_pair & representatives$new_representative=="YES" & representatives$study=="fromYimeng"),]

#final_matrix=p_matrix_depleted[ which(rownames(p_matrix_depleted) %in% family_pair_metadata$ID), ]

#remove uninteresting cells
#final_matrix=final_matrix[,-which( colSums((abs(final_matrix) < 50)) == nrow(final_matrix))]

#Spacing motifs: Convert the motif names (row labels) name + 5th column, Ns shortened
my_list=strsplit(rownames(final_matrix), "_")
# Find the length of the longest vector
max_length <- max(sapply(my_list, length))

# Pad shorter vectors with NA
padded_list <- lapply(my_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Convert the list to a data frame
# row_name_info <- do.call(rbind, padded_list)
# row_name_info[is.na(row_name_info[,9]),5]=row_name_info[is.na(row_name_info[,9]),6]
# row_name_info=as.data.frame(row_name_info[,-seq(6,10,1)])
# row_name_info$ID=rownames(final_matrix)

#Do this only for the spacing motifs
# Regular expression for a sequence of 'N's in the middle

#pattern <- "(?<!^)(NN+)(?!$)"

# Find positions of the pattern
# matches <- gregexpr(pattern, row_name_info[,5], perl = TRUE)
# matches_info=as.data.frame(cbind(row_name_info[,5],do.call(rbind, lapply(matches,attributes))))
# 
# replacement=paste0(unlist(matches_info$capture.length) ,"N")
# 
# matches_info$new_name=""

# for(i in 1:length(row_name_info[,5])){
#   matches_info$new_name[i]=gsub(pattern, replacement[i], row_name_info[i,5], perl = TRUE)
# }
# 
# row_name_info$new_name=matches_info$new_name
# 
# rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#Composite motifs
# row_name_info=family_pair_metadata
# 
# rownames(final_matrix)=paste0(family_pair_metadata$symbol, "_",family_pair_metadata$seed)


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


#fontsize=12


# bold=representatives$new_representative[match(row_name_info$ID, representatives$ID)]
# row_label_color=representatives$study[match(row_name_info$ID, representatives$ID)]
# row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")
# row_name_colors <- ifelse(row_label_color %in% c("fromYimeng"), "red", "black")

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

#Add the cell type group as annotation

result=as.data.frame(do.call(rbind, strsplit(family_pair_metadata$Lambert2018_families,"_")))
names(result)=c("X1", "X2")

#Cell type group names

#Colors for 15 cell type groups, use these in all figures

palette_colors <- brewer.pal(8, "Set2")

ct_group_colors=c("#66C2A5", #Stromal
                  "#FC8D62", #Islet
                  "#8DA0CB", #Brain glia
                  "#E78AC3", #Brain neuron
                  "#A6D854", #"Epithelial 1" 
                  "#FFD92F", #"Smooth Muscle" 
                  "#E5C494", #"Myocyte"        
                  "#B3B3B3", #  "Mural" 
                  "#E41A1C", #"Adrenal Cortical"
                  "#377EB8", #"Immune Lymp" 
                  "#FFFF33", #"Immune Myelo"  
                  "#A65628", #"Endothelial" 
                  "#000000", #"Epithelial 2" 
                  "#CCEBC5", #"Epithelial 3"
                  "#BEAED4") #"GI Epithelial"   
                  #"#FFFF99") 

names(ct_group_colors)=c( "Stromal", 
                          "Islet",
                          "Brain glia",
                          "Brain neuron", 
                          "Epithelial 1",
                          "Smooth Muscle",
                          "Myocyte",
                          "Mural",         
                          "Adrenal Cortical",
                          "Immune Lymp",
                          "Immune Myelo",
                          "Endothelial",
                          "Epithelial 2",
                          "Epithelial 3",
                          "GI Epithelial")          
           



# ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
#                                                     function(y) sapply(cell_type_groups, function(x) y %in% x)),
#                                              2, 
#                                              function(x) which(x==TRUE)) ]
# annotation_colors=ct_group_colors

#Sort according to ct_group names

#col_order=order(ct_group_names)
#final_matrix=final_matrix[,col_order]
#ct_group_names=ct_group_names[col_order]

# Generate a color palette
#num_colors <- length(unique(ct_group_names))
#palette_colors <- brewer.pal(num_colors, "Set2")  # Or any other palette

# Create a named vector of colors
#annotation_colors <- setNames(palette_colors, unique(ct_group_names))


# library(proxy)
# # Calculate cosine similarity
# cosine_similarity_rows <- proxy::simil(final_matrix, method = "cosine")
# # Convert to distance
# cosine_distance_rows <- 1 - cosine_similarity_rows
# 
# cosine_similarity_columns <- proxy::simil(t(final_matrix), method = "cosine")
# # Convert to distance
# cosine_distance_columns <- 1 - cosine_similarity_columns
# 
# dend1 = cluster_within_group(final_matrix, ct_group_names)


# pdf(file = paste0(scratch, "Figures/enrichment_",family_pair[1],".pdf"), width=14, 
#     height=12)
# 
# draw(ComplexHeatmap::Heatmap(final_matrix, 
#                         col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
#                         #col=colorRamp2(c(-100, 0, 100), c("navy", "white", "firebrick3")), 
#                         #cluster_rows = hclust(as.dist(cosine_distance_rows)), 
#                         #cluster_columns = FALSE, #hclust(as.dist(cosine_distance_columns)),
#                          cluster_rows = TRUE,
#                          cluster_columns = dend1, 
#                          clustering_distance_rows = "manhattan",
#                          clustering_method_rows = "complete",
#                          clustering_distance_columns = "pearson",
#                          clustering_method_columns = "complete",
#                         # 
#                         show_row_names = TRUE,
#                         show_column_names = TRUE,
#                         column_names_side = "top",
#                         row_names_side = "left",
#                         column_names_rot=45, #315,
#                       
#                         left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
#                                                       col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
#                                                       annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
#                                                                                    gp=gpar(fontsize=fontsize)),
#                                                       annotation_name_side="top",
#                                                       #annotation_name_rot=0,
#                                                       annotation_label=c("Protein family 1","Protein family 2"),
#                                                       #annotation_label="Protein\nfamily",
#                                                       gp=gpar(fontsize=fontsize),
#                                                       width = unit(10, "cm")),
#                         top_annotation=HeatmapAnnotation(`Cell type group`=ct_group_names,
#                                                          col=list(`Cell type group`=annotation_colors),
#                                                          annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
#                                                            gp=gpar(fontsize=fontsize)),
#                                                          annotation_name_side="right",
#                                                          #annotation_name_rot=0,
#                                                          #annotation_label=c("Protein family 1","Protein family 2"),
#                                                          #annotation_label="Protein\nfamily",
#                                                          gp=gpar(fontsize=fontsize)),
#                         show_column_dend = FALSE,
#                         show_row_dend = FALSE,
#                         row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
#                         height = nrow(final_matrix)*unit(0.618,"cm"), 
#                         width = ncol(final_matrix)*unit(1, "cm"),
#                         heatmap_legend_param=list(
#                           #title = "-Log_10(P-val)",
#                           #title = expression(-log[10](P-val)),
#                           title = expression(paste(-log[10],"(P-value)")),
#                           title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
#                           #title_position = "topleft",
#                           grid_height = unit(8, "mm"),
#                           grid_width = unit(8, "mm"),
#                           #tick_length = unit(0.8, "mm"),
#                           #border = NULL,
#                           #at = object@levels,
#                           #labels = at,
#                           at = c(0, 100),
#                           labels = c("0",">100 (Enriched)"),
#                           labels_gp = gpar(fontsize = fontsize),
#                           #labels_rot = 0,
#                           #nrow = NULL,
#                           #ncol = 1,
#                           #by_row = FALSE,
#                           legend_gp = gpar(fontsize = fontsize),
#                           #legend_height = NULL,
#                           #legend_width = NULL,
#                           legend_direction = "horizontal" #c("vertical", "horizontal"),
#                           # break_dist = NULL,
#                           #graphics = NULL,
#                           #param = NULL
#                         )
# ), heatmap_legend_side="bottom", annotation_legend_side = "bottom")
# 
# dev.off()


#Logistic regression coefficients

representatives$new_motif=FALSE

representatives$new_motif[representatives$study=="fromYimeng"]=TRUE 

 

regression_coefs=readRDS( paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/all_regression_coefs_maxscore.RDS"))

scaled_regression_coefs=scale(regression_coefs,center=TRUE, scale=TRUE)

final_matrix=regression_coefs[which(rownames(regression_coefs)%in% family_pair_metadata$ID),] #Only 4 representatives

#For spacing motifs
#rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#rownames(final_matrix)=row_name_info$new_name_full[match(rownames(final_matrix),row_name_info$ID)]

#Composites
ind=match(rownames(final_matrix),family_pair_metadata$ID)
family_pair_metadata=family_pair_metadata[ind,]
row_name_info=family_pair_metadata

rownames(final_matrix)=paste0(family_pair_metadata$symbol, "_",family_pair_metadata$seed)



final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 0.02)) ==nrow(final_matrix) )]



ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                     function(y) sapply(cell_type_groups, function(x) y %in% x)),
                              2, 
                              function(x) which(x==TRUE)) ]
#Sort according to ct_group names

col_order=order(ct_group_names)
final_matrix=final_matrix[,col_order]
ct_group_names=ct_group_names[col_order]

# Generate a color palette
#num_colors <- length(unique(ct_group_names))
#palette_colors <- brewer.pal(num_colors, "Set2")  # Or any other palette

# Create a named vector of colors
#annotation_colors <- setNames(palette_colors, unique(ct_group_names))

annotation_colors=ct_group_colors

#bold=representatives$new_representative[match(row_name_info$ID, representatives$ID)]
#row_label_color=representatives$study[match(row_name_info$ID, representatives$ID)]
#row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")
#row_name_colors <- ifelse(row_label_color %in% c("fromYimeng"), "red", "black")

result=as.data.frame(do.call(rbind, strsplit(family_pair_metadata$Lambert2018_families,"_")))
names(result)=c("X1", "X2")

fontsize=12

pdf(file = paste0(scratch, "Figures/regression_coeffs_",family_pair[1],".pdf"), width=16, 
    height=7)

draw(ComplexHeatmap::Heatmap(final_matrix, 
                             #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                             col=colorRamp2(c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), c("navy", "white", "firebrick3")), 
                             #col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")), 
                             cluster_rows = FALSE,
                             cluster_columns = FALSE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             column_names_rot=45, #315,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             clustering_distance_rows = "euclidean",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "euclidean",
                             clustering_method_columns = "complete",
                              left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
                                                            col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
                                                            annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                              gp=gpar(fontsize=fontsize)),
                                                            annotation_name_side="top",
                                                            #annotation_name_rot=0,
                                                            annotation_label=c("Protein family 1","Protein family 2"),
                                                            #annotation_label="Protein\nfamily",
                                                            gp=gpar(fontsize=fontsize),
                                                            width = unit(10, "cm")),
                             top_annotation=columnAnnotation(`Cell type group`=ct_group_names,
                                                             col=list(`Cell type group`=annotation_colors),
                                                             annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                             gp=gpar(fontsize=fontsize)),
                                                             annotation_name_side="right",
                                                             #annotation_name_rot=0,
                                                             #annotation_label=c("Protein family 1","Protein family 2"),
                                                             #annotation_label="Protein\nfamily",
                                                             gp=gpar(fontsize=fontsize)),
                                                             #width = unit(10, "cm")), 
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             height = nrow(final_matrix)*unit(0.618,"cm"), 
                             width = ncol(final_matrix)*unit(1, "cm"),
                             heatmap_legend_param=list(
                               #title = "-Log_10(P-val)",
                               #title = expression(-log[10](P-val)),
                               title = "Regression coefficient",
                               title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                               #title_position = "topleft",
                               grid_height = unit(8, "mm"),
                               grid_width = unit(8, "mm"),
                               #tick_length = unit(0.8, "mm"),
                               #border = NULL,
                               #at = object@levels,
                               #labels = at,
                               at = c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), #c(-0.05, 0, 0.05),
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
), heatmap_legend_side="bottom",annotation_legend_side = "bottom")

dev.off()


#Scaled regression coefficients

final_matrix=scaled_regression_coefs[which(rownames(scaled_regression_coefs)%in% family_pair_metadata$ID),]
final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 2)) ==nrow(final_matrix) )]

#Spacing motifs
#rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#Composites
ind=match(rownames(final_matrix),family_pair_metadata$ID)
family_pair_metadata=family_pair_metadata[ind,]
row_name_info=family_pair_metadata

rownames(final_matrix)=paste0(family_pair_metadata$symbol, "_",family_pair_metadata$seed)

ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]
#Sort according to ct_group names

col_order=order(ct_group_names)
final_matrix=final_matrix[,col_order]
ct_group_names=ct_group_names[col_order]

#Ets_Runt order

c("Endothelial", "Endothelial", "Endothelial", "Immune Lymp", "Immune Lymp","Immune Lymp", "Immune Lymp", "Immune Myelo",  "Immune Myelo", "Immune Myelo", 
"Adrenal Cortical", "Brain glia", "Epithelial 1", "Stromal") 

ind=c(3,4,5,7,8,9,10,11,12,13,1,2,6,14)
ct_group_names=ct_group_names[ind]
final_matrix=final_matrix[,ind]


#Forkhead_Ets order
ind=c(3,4,5,6,7,8,13,14,15,16,9,10,11,12,1,2,22,23,17,18,19)
ct_group_names=ct_group_names[ind]
final_matrix=final_matrix[,ind]

#Move goblet cells close to each other
ind_order=seq(1,ncol(final_matrix),1)
ind_order[13]=12
ind_order[12]=13
final_matrix=final_matrix[,ind_order] 

# Generate a color palette
#num_colors <- length(unique(ct_group_names))
#palette_colors <- brewer.pal(num_colors, "Set3")  # Or any other palette

# Create a named vector of colors
#annotation_colors <- setNames(palette_colors, unique(ct_group_names))

#row_label_color=representatives$study[match(row_name_info$ID, representatives$ID)]
#row_name_fontfaces <- ifelse(bold %in% c("YES"), "bold", "plain")
#row_name_colors <- ifelse(row_label_color %in% c("fromYimeng"), "red", "black")

#Do not draw protein family and seed

#rownames(final_matrix)=paste0(row_name_info$V1,"_",row_name_info$V2)
rownames(final_matrix)=row_name_info$symbol

#add motif

#Draw also reverse complements
image_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2_Logos/svg/", row_name_info$ID, ".svg")

# we only draw the image annotation for PNG images, while the others are the same
#install.packages("grImport2")
#install.packages("rsvg")

ha = HeatmapAnnotation(Logo = anno_image(image_svg, height= unit(8, "cm"), 
                                         border = FALSE,
                                         width= unit(8, "cm"),
                                         space = unit(0, "mm")), which="row") #, height=unit(0.618,"cm")) #, 


pdf(file = paste0(scratch, "Figures_version2.2/scaled_regression_coeffs_",family_pair[1],".pdf"), width=16, 
    height=7)



ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                             #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                             col=colorRamp2(c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), c("navy", "white", "firebrick3")), 
                             #col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")), 
                             cluster_rows = FALSE,
                             cluster_columns = FALSE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             column_names_rot=45, #315,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             clustering_distance_rows = "euclidean",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "euclidean",
                             clustering_method_columns = "complete",
                             # left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
                             #                               col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
                             #                               annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                             #                                 gp=gpar(fontsize=fontsize)),
                             #                               annotation_name_side="top",
                             #                               #annotation_name_rot=0,
                             #                               annotation_label=c("Protein family 1","Protein family 2"),
                             #                               #annotation_label="Protein\nfamily",
                             #                               gp=gpar(fontsize=fontsize),
                             #                               width = unit(10, "cm")),
                             top_annotation=columnAnnotation(`Cell type group`=ct_group_names,
                                                             col=list(`Cell type group`=annotation_colors),
                                                             annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                               gp=gpar(fontsize=fontsize)),
                                                             annotation_name_side="left",
                                                             #annotation_name_rot=0,
                                                             #annotation_label=c("Protein family 1","Protein family 2"),
                                                             #annotation_label="Protein\nfamily",
                                                             gp=gpar(fontsize=fontsize)),
                             #width = unit(10, "cm")), 
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             height = nrow(final_matrix)*unit(1,"cm"), 
                             width = ncol(final_matrix)*unit(0.8, "cm"),
                             heatmap_legend_param=list(
                               #title = "-Log_10(P-val)",
                               #title = expression(-log[10](P-val)),
                               title = "Scaled regression coefficient",
                               title_gp = gpar(fontsize = fontsize),
                               #title_position = "topleft",
                               grid_height = unit(8, "mm"),
                               grid_width = unit(8, "mm"),
                               #tick_length = unit(0.8, "mm"),
                               #border = NULL,
                               #at = object@levels,
                               #labels = at,
                               at = c(-floor(max(abs(final_matrix))) , 0, floor(max(abs(final_matrix))) ),
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



draw(ht_list+ha, 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.off()


#Draw both family pairs together


family_pair=c("Ets_Runt", "Runt_Ets", "Forkhead_Ets", "Ets_ForkHead") 

family_pair_metadata=representatives[which(representatives$Lambert2018_families %in% family_pair & representatives$new_representative=="YES" & representatives$study=="fromYimeng"),] #7

final_matrix=regression_coefs[which(rownames(regression_coefs)%in% family_pair_metadata$ID),] 


#Composites
ind=match(rownames(final_matrix),family_pair_metadata$ID)
family_pair_metadata=family_pair_metadata[ind,]
row_name_info=family_pair_metadata

rownames(final_matrix)=paste0(family_pair_metadata$symbol, "_",family_pair_metadata$seed)


#remove uninteresting cells
#final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 0.02)) ==nrow(final_matrix) )]
final_matrix=final_matrix[,-which( colSums( ( abs(final_matrix) < 0.02| (final_matrix < 0))) ==nrow(final_matrix) )]




ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]
#Sort according to ct_group names

col_order=order(ct_group_names)
final_matrix=final_matrix[,col_order]
ct_group_names=ct_group_names[col_order]

annotation_colors=ct_group_colors
result=as.data.frame(do.call(rbind, strsplit(family_pair_metadata$Lambert2018_families,"_")))
names(result)=c("X1", "X2")

fontsize=12

pdf(file = paste0(scratch, "Figures/regression_coeffs_",family_pair[1],".pdf"), width=16, 
    height=7)

draw(ComplexHeatmap::Heatmap(final_matrix, 
                             #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                             col=colorRamp2(c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), c("navy", "white", "firebrick3")), 
                             #col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")), 
                             cluster_rows = FALSE,
                             cluster_columns = FALSE, 
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             column_names_side = "top",
                             row_names_side = "left",
                             column_names_rot=45, #315,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             clustering_distance_rows = "pearson",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "pearson",
                             clustering_method_columns = "complete",
                             left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
                                                           col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
                                                           annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                             gp=gpar(fontsize=fontsize)),
                                                           annotation_name_side="top",
                                                           #annotation_name_rot=0,
                                                           annotation_label=c("Protein family 1","Protein family 2"),
                                                           #annotation_label="Protein\nfamily",
                                                           gp=gpar(fontsize=fontsize),
                                                           width = unit(10, "cm")),
                             top_annotation=columnAnnotation(`Cell type group`=ct_group_names,
                                                             col=list(`Cell type group`=annotation_colors),
                                                             annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                               gp=gpar(fontsize=fontsize)),
                                                             annotation_name_side="right",
                                                             #annotation_name_rot=0,
                                                             #annotation_label=c("Protein family 1","Protein family 2"),
                                                             #annotation_label="Protein\nfamily",
                                                             gp=gpar(fontsize=fontsize)),
                             #width = unit(10, "cm")), 
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
                             #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                             height = nrow(final_matrix)*unit(0.618,"cm"), 
                             width = ncol(final_matrix)*unit(1, "cm"),
                             heatmap_legend_param=list(
                               #title = "-Log_10(P-val)",
                               #title = expression(-log[10](P-val)),
                               title = "Regression coefficient",
                               title_gp = gpar(fontsize = fontsize+2, fontface = "bold"),
                               #title_position = "topleft",
                               grid_height = unit(8, "mm"),
                               grid_width = unit(8, "mm"),
                               #tick_length = unit(0.8, "mm"),
                               #border = NULL,
                               #at = object@levels,
                               #labels = at,
                               at = c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), #c(-0.05, 0, 0.05),
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
), heatmap_legend_side="bottom",annotation_legend_side = "bottom")

dev.off()


#Scaled regression coefficients

final_matrix=scaled_regression_coefs[which(rownames(scaled_regression_coefs)%in% family_pair_metadata$ID),]
#final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 2)) ==nrow(final_matrix) )]

final_matrix=final_matrix[,-which( colSums( (abs(final_matrix) < 2| (final_matrix < 0))) ==nrow(final_matrix) )]


#Spacing motifs
#rownames(final_matrix)=apply(row_name_info[,c(1,2,7)],1,paste0,collapse="_")

#Composites
ind=match(rownames(final_matrix),family_pair_metadata$ID)
family_pair_metadata=family_pair_metadata[ind,]
row_name_info=family_pair_metadata

rownames(final_matrix)=paste0(family_pair_metadata$symbol, "_",family_pair_metadata$seed)

ct_group_names=names(cell_type_groups)[apply(sapply(colnames(final_matrix), 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]
#Sort according to ct_group names

col_order=order(ct_group_names)
final_matrix=final_matrix[,col_order]
ct_group_names=ct_group_names[col_order]

#Ets_Runt_Forkhead order

ind=c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,1,2,17,18)
ct_group_names=ct_group_names[ind]
final_matrix=final_matrix[,ind]


#rownames(final_matrix)=paste0(row_name_info$V1,"_",row_name_info$V2)
rownames(final_matrix)=row_name_info$symbol

#add motif

#Draw also reverse complements
image_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2_Logos/svg/", row_name_info$ID, ".svg")

# we only draw the image annotation for PNG images, while the others are the same
#install.packages("grImport2")
#install.packages("rsvg")

ha = HeatmapAnnotation(Logo = anno_image(image_svg, height= unit(8, "cm"), 
                                         border = FALSE,
                                         width= unit(8, "cm"),
                                         space = unit(0, "mm")), which="row") #, height=unit(0.618,"cm")) #, 


pdf(file = paste0(scratch, "Figures_version2.2/scaled_regression_coeffs_Ets_Runt_Forkhead_positive.pdf"), width=16, 
    height=10)



ht_list=ComplexHeatmap::Heatmap(final_matrix, 
                                #col=colorRamp2(c( 0, 100), c("white", "firebrick3")), 
                                col=colorRamp2(c(-ceiling(max(abs(final_matrix))*100)/100 , 0, ceiling(max(abs(final_matrix))*100)/100 ), c("navy", "white", "firebrick3")), 
                                #col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")), 
                                cluster_rows = FALSE,
                                cluster_columns = FALSE, 
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                column_names_side = "top",
                                row_names_side = "left",
                                column_names_rot=45, #315,
                                #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                                clustering_distance_rows = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_distance_columns = "euclidean",
                                clustering_method_columns = "complete",
                                left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
                                                              col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
                                                              annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                                gp=gpar(fontsize=fontsize)),
                                                              annotation_name_side="top",
                                                              #annotation_name_rot=0,
                                                              annotation_label=c("Protein family 1","Protein family 2"),
                                                              #annotation_label="Protein\nfamily",
                                                              gp=gpar(fontsize=fontsize),
                                                              width = unit(10, "cm")),
                                top_annotation=columnAnnotation(`Cell-type group`=ct_group_names,
                                                                col=list(`Cell-type group`=annotation_colors),
                                                                annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
                                                                  gp=gpar(fontsize=fontsize)),
                                                                annotation_name_side="right",
                                                                #annotation_name_rot=0,
                                                                #annotation_label=c("Protein family 1","Protein family 2"),
                                                                #annotation_label="Protein\nfamily",
                                                                gp=gpar(fontsize=fontsize)),
                                #width = unit(10, "cm")), 
                                show_column_dend = FALSE,
                                show_row_dend = FALSE,
                                #row_names_gp = gpar(col = row_name_colors, fontface=row_name_fontfaces),
                                height = nrow(final_matrix)*unit(1,"cm"), 
                                width = ncol(final_matrix)*unit(0.8, "cm"),
                                heatmap_legend_param=list(
                                  #title = "-Log_10(P-val)",
                                  #title = expression(-log[10](P-val)),
                                  title = "Scaled regression coefficient",
                                  title_gp = gpar(fontsize = fontsize),
                                  #title_position = "topleft",
                                  grid_height = unit(8, "mm"),
                                  grid_width = unit(8, "mm"),
                                  #tick_length = unit(0.8, "mm"),
                                  #border = NULL,
                                  #at = object@levels,
                                  #labels = at,
                                  at = c(floor(min(final_matrix)) , 0, ceiling(max(abs(final_matrix))) ),
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



draw(ht_list+ha, 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.off()

