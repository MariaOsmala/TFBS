library("rtracklayer")
library("dbplyr")
library("dplyr")
library("tidyverse")
library("pheatmap")
library("cluster")

#short names for TFs


scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

data_path="/scratch/project_2006203/"

representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives_IC_length_better_ID_filenames_match_numbers_MOODS_threshold.tsv", sep="\t", header=TRUE)

#test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
test=representatives$ID
rep_motifs=test[which(representatives$new_representative=="YES")]

#These are log_10(p-values)
p_matrix=readRDS(  file = paste0(scratch, "RData/Teemu_p_matrix_human_cell_type_restricted.Rds")) #

shorter_names=representatives %>% filter(test %in% rownames(p_matrix)) %>% select(symbol, ID, experiment, new_representative, Lambert2018.families)

#These are log2( (k / n) / (m / N) )
e_matrix=readRDS(  file = paste0(scratch, "RData/Teemu_e_matrix_human_cell_type_restricted.Rds")) #

#Heatmaps showing Gene Ontology TF motifs with maximal enrichment in cell-type-restricted cCREs of selected cell types. 
#Only the most enriched TF motif in each of the previously identified motif archetypes (Vierstra et al., 2020) was selected as the representative and 
#the top 10 motifs were selected for each cell type. Color represents −log10P.

#which are enriched 
infinites=which(is.infinite(e_matrix))
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
p_matrix_depleted=-p_matrix
p_matrix_depleted[depleted]=-p_matrix_depleted[depleted]

#plot heatmap



min(p_matrix_depleted) #-264.5894
max(p_matrix_depleted) #1578.851



#Epithelial
#Visualise the Hepatocyte, Acinar, Parietal, Chief, Foveolar, G. Neuroendo, Ductal, 

#which are depleted 

p_matrix_depleted_orig=p_matrix_depleted


rownames(p_matrix_depleted)=shorter_names$symbol
rownames(e_matrix)=shorter_names$symbol

#Different cell-lineage groups



epithelial=p_matrix_depleted[, which( colnames(p_matrix_depleted) %in% 
                                        c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal") )]
e_epithelial=e_matrix[, which( colnames(e_matrix) %in% 
                                          c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal") )]

epithelial_orig=p_matrix_depleted_orig[, which( colnames(p_matrix_depleted_orig) %in% 
                                        c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal") )]


       
#what are the 10 most enriched motifs for these cell lines       

most_enriched<-list()
most_enriched_orig<-list()
for(c in colnames(epithelial)){
  
  most_enriched[[c]]=names(sort(epithelial[,c], decreasing = TRUE)[1:10])
  most_enriched_orig[[c]]=names(sort(epithelial_orig[,c], decreasing = TRUE)[1:10])
  
}

length(unlist(most_enriched)) #70
most_enriched=unique(unlist(most_enriched)) #33

length(unlist(most_enriched_orig)) #70
most_enriched_orig=unique(unlist(most_enriched_orig)) #44

#index
index=which(rownames(epithelial_orig) %in% unlist(most_enriched_orig))

epithelial=epithelial[index,] 
e_epithelial=e_epithelial[index,] 
epithelial_orig=epithelial_orig[index,] 


#epithelial[
#  order( epithelial[,1], epithelial[,2],epithelial[,3],epithelial[,4],epithelial[,5],epithelial[,6] , decreasing=FALSE),
#]

#epithelial[
#  with(epithelial, order("Foveolar", decreasing=FALSE)),
#]

#epithelial=epithelial[
#  order( epithelial[,c( "Foveolar")] ),
#]

epithelial=epithelial[,c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal")]
epithelial_orig=epithelial_orig[,c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal")]
e_epithelial=e_epithelial[,c("Hepatocyte", "Acinar", "Parietal", "Chief", "Foveolar", "Gastric Neuroendocrine", "Ductal")]

#epithelial=epithelial[order(rownames(epithelial)),]

#Vierstra
#epithelial=epithelial[c("HNF1A","HNF4A","ESRRA","GATA2","GATA3","GATA5", "FOXA1", "FOXA2" , "FOXC2",  "FOXD2",  "ONECUT1",  "ONECUT2",  "NR2F1", "BARHL1_TEAD4", "CEBPB_FOXD2",
#  "CEBPG",  "FIGLA",  "FOSL1", "GSC2_TEAD4", "HOXA10_TEAD4", "HOXB2_TCF3", "HOXD12_TEAD4", "ID4", "JUN", "LMX1A_BACH2",  "MYBL1_FIGLA",  "NEUROG1_TCF3",
#"NR2C1", "OTX1_TEAD4", "POU4F1", "RARA","RFX1", "RFX3",  "RFX3_FIGLA", "TCF12", "TEAD2", "TFAP4_FLI1" ),]

#epithelial=epithelial[c("HNF1A", "HNF4A", "ESRRA", "GATA2", "GATA3", "GATA5",  "FOXA1",  "FOXA2",  "FOXD2",  "ONECUT1", "ONECUT2",  "CEBPB_FOXD2", "CEBPG", "FIGLA",  "GSC2_TEAD4", 
#                        "HOXB2_TCF3",  "ID4", "JUN",  
#"LMX1A_BACH2",  "MYBL1_FIGLA",  "NEUROG1_TCF3", "NR2C1",  "POU4F1",   "RARA",  "RFX1",  "RFX3", "RFX3_FIGLA",  "TCF12",  "TEAD1",  "TEAD2",  "ETV5_TCF3", "FOXC2_TCF3", "CUX1", "CUX2", "SOX10", "KLF12", "RFX7"),]                 
                        
#e_epithelial=e_epithelial[c("HNF1A", "HNF4A", "ESRRA", "GATA2", "GATA3", "GATA5",  "FOXA1",  "FOXA2",  "FOXD2",  "ONECUT1", "ONECUT2",  "CEBPB_FOXD2", "CEBPG", "FIGLA",  "GSC2_TEAD4", 
#                        "HOXB2_TCF3",  "ID4", "JUN",  
#                        "LMX1A_BACH2",  "MYBL1_FIGLA",  "NEUROG1_TCF3", "NR2C1",  "POU4F1",   "RARA",  "RFX1",  "RFX3", "RFX3_FIGLA",  "TCF12",  "TEAD1",  "TEAD2",  "ETV5_TCF3", "FOXC2_TCF3",
#                        "CUX1", "CUX2", "SOX10", "KLF12", "RFX7"),]                 

tmp=c("HNF1A_HT-SELEX_TGAGCA20NCGA_W_NRTTAATNATTAACN_1_2_YES", 
"HNF4A_HT-SELEX_TACCTT40NCGA_KO_NRGTCCAAAGGTCRN_1_3_NO",
"HNF4A_HT-SELEX_TACCTT40NCGA_KO_NRGTCCAAAGTCCRN_1_3_NO" ,                        
"HNF4A_HT-SELEX_TCCGTG40NTGC_AI_RRGGTCAAAGTCCRNN_1_3_YES",
"HNF4A_Methyl-HT-SELEX_TTGTTT40NGGC_KO_NRGGTCAAAGGTCAN_1_3_NO"  ,
"ESRRA_HT-SELEX_TAGACC40NAGT_KO_NYCAAGGTCAN_1_3_NO",
"GATA2_HT-SELEX_TGAGGT40NTCC_KR_NYGATAASN_1_2_YES",
"GATA2_Methyl-HT-SELEX_TATGCA40NGAG_KR_NWGATAASN_1_2_YES",
"GATA3_HT-SELEX_TTCGTT40NGCC_KX_NGATAAGATCW_1_3_YES",
"GATA5_Methyl-HT-SELEX_TTACTT40NTAT_KW_NWGATAASRN_1_2_NO" ,
"FOXA1_HT-SELEX_TTCTAA40NAAT_KN_TRNGTAAACA_1_3b1_NA",                            
"FOXA2_HT-SELEX_TAGCTG40NCTA_KT_NWNWGTMAATATTKRYNYWN_2_4_NO",
"FOXC2_TCF3_TGACGT40NAGC_YIII_NCACCTGNRTAAAYAN_m1_c3b0u_short_composite",
"FOXD2_HT-SELEX_TCGCCG40NTGC_KO_NYWANGTAAACAN_1_4_YES",
"FOXD2_Methyl-HT-SELEX_TCATCT40NATT_KO_NYWANGTAAACAN_1_4_YES"   ,
"ONECUT1_Methyl-HT-SELEX_TGGCCG40NGCT_KW_NTATTGATCSGN_1_4_NO",                   
"ONECUT2_Methyl-HT-SELEX_TCCGGA40NGTC_KZ_NTATTGATTWN_1_2_NO" ,
"CEBPB_FOXD2_TCCTCT40NGTC_YWI_NTTRYGTAAACAN_m1_c3b0_short_composite",
"CEBPG_HT-SELEX_TAAAAT20NCG_AC_NTTRCGCAAY_1_2_NO",
"FIGLA_HT-SELEX_AGATA14N_U_NNCACCTGNN_1_4_NO" ,
"GSC2_TEAD4_TTTCTA40NCAT_YJII_NMRTTAACATTCCN_m1_c3b0_short_spacing" ,
"HOXB2_TCF3_CAP-SELEX_TGAAGC40NCTA_AY_NCACCTGNNNNNMATTA_1_3_YES" ,
"ID4_HT-SELEX_CTACA14N_U_NRCACCTGNN_1_4_YES"  ,
"JUN_HT-SELEX_TCAGTG40NGAT_KAE_NATGACKCATN_1_3b0_NO",
"LMX1A_BACH2_TTATTT40NCGT_YACIIII_TGASTCANNNNNNNNNNTAATTA_m1_c3b0_short_spacing",
"MYBL1_FIGLA_CAP-SELEX_TCAGCC40NTTC_AX_NCASSTGNNNNNNNCSGTTR_1_3_YES" ,
"NEUROG1_TCF3_TCGCAC40NCCT_YIII_NMCATCTGYYN_m1_c3b0u_short_composite",
"NR2C1_HT-SELEX_TATTCA40NAGT_KT_NRRGGTCAN_1_4_YES",                              
"NR2C1_Methyl-HT-SELEX_TACGGC40NAGT_KT_NRAGGTCAN_1_4_YES",
"POU4F1_HT-SELEX_TACAAA20NAAC_Y_ATGMATAATTAATG_2_3_YES" , 
"RARA_HT-SELEX_TGAACC40NCCA_KS_NRRGGTCANN_1_4_NO" ,  
"RFX1_HT-SELEX_TTTAGG40NTCT_KS_NGTTRCCATGGYAACN_1_3_NO" ,
"RFX3_HT-SELEX_TCCACC40NTTA_KS_NGTTGCCWAGCAACN_1_2_NO",
"RFX3_FIGLA_CAP-SELEX_TCGGCA40NAGC_AY_TRGYAACNNNNCASSTGNN_1_3_YES",
"RFX3_FIGLA_CAP-SELEX_TCGGCA40NAGC_AY_GTTGCYNNNNNNNNNNNNCASSTG_1_3_YES" ,
"RFX7_HT-SELEX_TGAGTT40NCGG_KT_NGTTGCYAN_1_3_NO",
"TCF12_HT-SELEX_TTTGAG40NATT_KS_NCACSTGN_1_3_NO",
"TEAD1_HT-SELEX_TCTTAG20NATG_W_RCATTCCNNRCATWCCN_2_4_YES" ,
"TEAD2_Methyl-HT-SELEX_TCGCAA40NTGT_KX_NRCATTCCWN_1_2_NO"  ,
"ETV5_TCF3_CAP-SELEX_TCGTCC40NGCC_AY_RNCGGAAGNNNNNCASSTGN_1_2_YES",
"CUX1_HT-SELEX_TTAGCG40NAGA_KP_NYATYGATYN_1_3_YES",        
"CUX2_HT-SELEX_TCCACC40NTTA_KP_NTGATCGATYRN_1_3_NO",  
"KLF12_Methyl-HT-SELEX_TAACTG40NGTA_KX_NRCCACGCCCW_1_3_YES",  "SOX10_HT-SELEX_TTGCCA40NGTG_KAF_AACAATNNNNNNATTGTT_1_3_YES")     

                  

#rownames(epithelial_orig)[match(tmp,rownames(epithelial_orig))]
epithelial=epithelial[match(tmp,rownames(epithelial_orig)),]
e_epithelial=e_epithelial[match(tmp,rownames(epithelial_orig)),]                      
epithelial_orig=epithelial_orig[match(tmp,rownames(epithelial_orig)),]                      


#save.image(file= "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/epithelial.RData")
load(file= "/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/epithelial.RData")


library(grid)
source("/scratch/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R")

library("RColorBrewer")
paletteLength=100
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdBu")))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))

library("pheatmap")
pdf(file = paste0(scratch, "Figures/top300000/adult-epithelial-heatmap.pdf"), width=8, height=4 )

pheatmap(t(epithelial), color=color, clustering_method="complete",legend=TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks= myBreaks, row_names_side = "left", row_dend_side = "right")
dev.off()

#Interactive heatmap

library(heatmaply)
library(RColorBrewer)
#install.packages("heatmaply")
epithelial_saturated=epithelial
saturated_threshold=100
epithelial_saturated[epithelial_saturated > saturated_threshold] <- saturated_threshold
epithelial_saturated[epithelial_saturated < -saturated_threshold] <- -saturated_threshold

#get the info for the TFs from metadata table

metadata=representatives[match(rownames(epithelial_orig), representatives$ID),]
[#1] "ID"                   "symbol"               "clone"                "family"               "Lambert2018.families"
#[6] "organism"             "study"                "experiment"           "ligand"               "batch"               
#[11] "seed"                 "multinomial"          "cycle"                "representative"       "new_representative"  
#[16] "short"                "type"                 "comment"              "filename"             "IC"                  
#[21] "length" "match_numbers", "MOODS_threshold"
  
  test=matrix(rep(metadata$ID,each=7), ncol=44, byrow=FALSE)

mat=t(epithelial_orig) #7 x 44
mat[]<-paste("Motif: ", matrix(rep(metadata$ID,7), ncol = 44, byrow=TRUE) ,
             "\nFamily: ", matrix(rep(metadata$Lambert2018.families,7), ncol = 44, byrow=TRUE),
             "\nStudy: ", matrix(rep(metadata$study,7), ncol = 44, byrow=TRUE),
             "\n-log10 pvalue:", mat, 
             "\nlog2FC: ", t(e_epithelial),
             "\nIC: ", matrix(rep(metadata$IC,7), ncol = 44, byrow=TRUE),
             "\nlength: ", matrix(rep(metadata$length,7), ncol = 44, byrow=TRUE),
             "\nmatch numbers: ", matrix(rep(metadata$match_numbers,7), ncol = 44, byrow=TRUE),
             "\nMOODS threshold: ", matrix(rep(metadata$MOODS_threshold,7), ncol = 44, byrow=TRUE)
             )

#mat[] <- 
#  
#  lapply( colnames(mat), function(colname) {
#  paste0("\nFamily ", metadata[ metadata$ID==colname, "Lambert2018.families"])} )
# Create heatmap with heatmaply
heatmaply(t(epithelial_saturated),
          col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
          legend=TRUE,
          showticklabels = T,
          ColSideColors = data.frame(Family=metadata$Lambert2018.families),
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
          main = "Epithelial cells",  # Title of the plot
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
          file="/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys/epithelial.html"
)



browseURL("/scratch/project_2006203/TFBS/ATAC-seq-peaks/heatmaplys/epithelial.html")










library("RColorBrewer")
paletteLength=100
color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdBu")))(paletteLength)

#color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-100, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(100/paletteLength, max(100), length.out=floor(paletteLength/2)))


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




library(heatmaply)
library(RColorBrewer)

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





