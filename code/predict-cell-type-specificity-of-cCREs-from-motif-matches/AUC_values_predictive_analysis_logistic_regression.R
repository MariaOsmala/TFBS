library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")

#install.packages("ggbreak")
#library("GRanges")
library("readr")
library("caret")
library("glmnet")
library("reshape2")
library(ggrepel) #
library("gridExtra")
library("openxlsx")



.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")
library(ggbreak)


#The logistic regression model predicts whether a cCRE i is specific to certain cell type (Class 1) or not (0). 
#The model is trained for all 111 cell types.

#In the simplest model, the covariates are presence/absence (0 or 1)  or counts (0,1,2,...)  
#of motif matches (individual representative PWMS ~1000 features). 
#Or use the binding affinities of motif matches? Positional logistic regression classifier? 
#Use R-packages caret and glmnet. 


setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
data_path="/scratch/project_2006203/TFBS/"

scratch="/scratch/project_2006203/TFBS/ATAC-seq-peaks/"

representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

representatives$new_motif=FALSE
representatives$new_motif[representatives$study=="fromYimeng"]=TRUE   
length(unique(rep_motifs))
#1232

#representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")

#motifs=names(representative_motif_matches)
#NOT NEEDED HERE
#cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
#length(cCREs_list)

#sum(sapply(cCREs_list, length)) #1 091 805

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"

#names(cCREs_list) %in% Adult_Celltypes_mapping$celltype
#Adult_Celltypes_mapping$celltype %in% names(cCREs_list) 

#test=data.frame( names(cCREs_list), Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ])

#names(cCREs_list)=Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ]
#all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
#all_cCREs=unique(all_cCREs)
#N=length(all_cCREs) #all possible cCREs #435142


#Cell type groups from Zhang et al. Figure 1

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



# aucs_ctg=list()
# assessments_ctg=list()
# cvfits_ctg=list()
# cvfit_final_ctg=list()
# 
# aucs_maxscore_ctg=list()
# assessments_maxscore_ctg=list()
# cvfits_maxscore_ctg=list()
# cvfit_maxscore_final_ctg=list()

aucs_test_all=data.frame()

for( cell_type_group in names(cell_type_groups) ){ 

  #cell_type_group=names(cell_type_groups)[1]
  print(cell_type_group)
  load(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/logistic_regression_processed_",cell_type_group,".RData"))
  
  aucs_test$`Celltype group`=cell_type_group
  
  aucs_test = aucs_test %>% select(id, `Celltype group`, everything())

  aucs_test_all=rbind(aucs_test_all, aucs_test)
    
}
  
saveRDS(aucs_test_all, paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/AUCs_all.RData"))

#Compute mean and standard deviation for each Cell type, need to add cell type group info
aucs_summary_cell_types <- aucs_test_all %>%
  group_by(`Cell type`, method) %>%
  summarize(
    mean_value = mean(AUC, na.rm = TRUE),
    sd_value = sd(AUC, na.rm = TRUE)
  )

aucs_summary_cell_types$`Celltype group`=names(cell_type_groups)[apply(sapply(aucs_summary_cell_types$`Cell type`, 
                                                    function(y) sapply(cell_type_groups, function(x) y %in% x)),
                                             2, 
                                             function(x) which(x==TRUE)) ]
aucs_summary_cell_types = aucs_summary_cell_types %>% select(`Celltype group`, everything())

#Compute mean and standard deviation for each Cell type group, order the results based on this (maxscore)
aucs_summary_cell_type_groups <- aucs_test_all %>%
  group_by(`Celltype group`, method) %>%
  summarize(
    mean_value = mean(AUC, na.rm = TRUE),
    sd_value = sd(AUC, na.rm = TRUE)
  ) %>% filter(method=="score") %>% arrange(desc(mean_value))


# Assuming your data frame is named my_data_frame
# and the column you want to convert is named my_column




ct_group_order=unique(aucs_summary_cell_type_groups$`Celltype group`)
aucs_summary_cell_types$`Celltype group` <- factor(aucs_summary_cell_types$`Celltype group`, levels = ct_group_order)

aucs_summary_cell_types <- aucs_summary_cell_types[order(aucs_summary_cell_types$`Celltype group`), ]

library(tidyverse)
aucs_summary_cell_types=aucs_summary_cell_types %>% pivot_wider(names_from=method, values_from=c(mean_value, sd_value))

aucs_summary_cell_types = aucs_summary_cell_types %>% select(`Celltype group`, `Cell type`, mean_value_score, sd_value_score, mean_value_presence, sd_value_presence)

#Write table in excel and highlight the higher mean
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
# Assuming your data frame is named my_data_frame
writeData(wb, "Sheet1", aucs_summary_cell_types)
# Assume column A (index 1) and column B (index 2) are the columns to compare
# Loop through the rows of the data frame
for (row in 1:nrow(aucs_summary_cell_types)) {
  if (aucs_summary_cell_types[row, "mean_value_score"] > aucs_summary_cell_types[row, "mean_value_presence"]) {
    # Apply bold formatting to the cell in column A for this row
    addStyle(wb, "Sheet1", style = createStyle(textDecoration = "bold"), rows = row + 1, cols = 3, gridExpand = FALSE)
  }else{
    addStyle(wb, "Sheet1", style = createStyle(textDecoration = "bold"), rows = row + 1, cols = 5, gridExpand = FALSE)
  }
}

saveWorkbook(wb, "AUCS_version2.2.xlsx", overwrite = TRUE)


library("broom")
#t-test for the difference between the maxscore and presence (highlight the significantly higher, or just highlight the mean)
test=aucs_test_all  %>% pivot_wider(names_from=method, values_from=AUC)

result <-test  %>%
  group_by(`Cell type`) %>%
  do(tidy(t.test(.$presence, .$score, paired=F, alternative="two.sided") ))

result %>% filter(p.value < 0.05)

result <-test  %>%
  group_by(`Cell type`) %>%
  do(tidy(t.test( .$score, .$presence, paired=F, alternative="greater") ))

aucs_summary_cell_types %>% filter(`Cell type` %in% c("Basal Epithelial (Mammary)",      
                                   "Foveolar Cell",
                                   "Fibroblast (Peripheral Nerve" ))


t.test(test$presence,test$score , paired=F, alternative = "two.sided")


# In the second step, we run another t-test as:
#   
#   (stutest <- t.test( all_auc, binding_auc, paired=T, alternative = c("less")))
# 
# Paired t-test
# 
# data:  all_auc and binding_auc
# t = -11.555, df = 99, p-value < 2.2e-16
# alternative hypothesis: true difference in means is less than 0
# 95 percent confidence interval:
#   -Inf -0.006235253
# sample estimates:
#   mean of the differences 
# -0.007281619 
# Here, since p-value is less than 0.05, significance level, we can reject H0 and accept H1. 
#Here H1 says that mean of samples in the first group is less than mean of samples in the second group.
# 


pdf(file = paste0(scratch, "Figures/cell_group_AUC_boxplots/",cell_type_group,".pdf"), width=12*(length(cell_type_groups[[cell_type_group]])/7), height=4 )
  
gg=ggplot(aucs_test, aes(x=factor(`Cell type`), y=AUC, fill=method)) + 
      geom_boxplot() + 
      labs(title="Cell-type-specific cCRE prediction based on motif matches, cross-validation,\nPerformance of logistic regression on test data",
           x="Cell type ",y="AUC") +
      theme_minimal()+labs(title=cell_type_group)
  plot(gg)
  
  dev.off()
  
  #write excel
  
  
  # Create a new workbook
  wb = createWorkbook()
  
  # Add a worksheet
  addWorksheet(wb, "presence")

  #Importance of the regression weights
  #Presence
  plots=list()
  plots_new=list()
  options(ggrepel.max.overlaps = Inf)
  
  i=1
  
  for(cell_type in names(cvfit_final_list) ){
    #cell_type=names(cvfit_final_list)[1]  
    regression_coef=data.frame(name=rownames(coef(cvfit_final_list[[cell_type]], s = "lambda.1se"))[-1], 
                               coef=as.matrix(coef(cvfit_final_list[[cell_type]], s = "lambda.1se"))[-1] )
    
    # 
    regression_coef=regression_coef[order(regression_coef$coef, decreasing =TRUE),]
    # 
    regression_coef$order=1:nrow(regression_coef)
    # 
    regression_coef$symbol=representatives$symbol[match(regression_coef$name, representatives$ID)]
    
    regression_coef$new_motif=representatives$new_motif[match(regression_coef$name, representatives$ID)]
    
    # 
    # #Plot regression coefficient as the function of the feature importance
    
    writeData(wb, "presence", x=cell_type, startCol=(i-1)*5+1, startRow=1)
    addStyle(wb, "presence", style=createStyle(border="left"), cols=(i-1)*5+1, rows=1)
    writeData(wb, "presence", x=cell_type, startCol=(i-1)*5+2, startRow=1)
    writeData(wb, "presence", x=cell_type, startCol=(i-1)*5+3, startRow=1)
    writeData(wb, "presence", x=cell_type, startCol=(i-1)*5+4, startRow=1)
    writeData(wb, "presence", x=cell_type, startCol=(i-1)*5+5, startRow=1)
    addStyle(wb, "presence", style=createStyle(border="right"), cols=(i-1)*5+5, rows=1)
    
    writeData(wb, "presence", x=regression_coef, startCol=(i-1)*5+1, startRow=2, 
              borders="surrounding", headerStyle = createStyle(border=c("left","right")) )
    
    # Define the style for the rows you want to color (e.g., red background)
    style_red <- createStyle(fgFill = "red")  # FF0000 is the hex code for red
    
    # Apply the style to the desired rows (e.g., rows 2 and 4)
    addStyle(wb, "presence", style = style_red, rows = 2+which(regression_coef$new_motif==TRUE), cols = ((i-1)*5+1):((i-1)*5+5), gridExpand = TRUE)
    
   
    
    # 
    plots[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point(aes(color = new_motif))+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      #geom_text_repel(data=subset(significant_motifs, new_motif), aes(label=shorter_name ) ,color="red"
      geom_text_repel(data=subset(regression_coef, (coef > 0.25 | coef < -0.25) & !new_motif, max.overlaps=Inf), 
                                                                        aes(label=symbol),color="black", box.padding=0.5)+
      geom_text_repel(data=subset(regression_coef, (coef > 0.25 | coef < -0.25) & new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="red", box.padding=0.5)+
      #geom_text_repel(data=subset(regression_coef, coef < -0.25, max.overlaps=Inf), 
      #                aes(label=symbol), box.padding=0.5)+
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient", color="New motif")+
      theme_minimal()
    
    plots_new[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point(aes(color = new_motif))+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      #geom_text_repel(data=subset(significant_motifs, new_motif), aes(label=shorter_name ) ,color="red"
      #geom_text_repel(data=subset(regression_coef, (coef > 0.25 | coef < -0.25) & !new_motif, max.overlaps=Inf), 
      #                aes(label=symbol),color="black", box.padding=0.6)+
      geom_text_repel(data=subset(regression_coef, (coef > 0.1 ) & new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="red", box.padding=unit(0.5, "lines"))+
      #geom_text_repel(data=subset(regression_coef, coef < -0.25, max.overlaps=Inf), 
      #                aes(label=symbol), box.padding=0.5)+
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient", color="New motif")+
      theme_minimal()
    
    
  # 
    i=i+1
  }
  
  
  
 
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list))))
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs_new_motifs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  
  # Add a worksheet
  addWorksheet(wb, "scores")
  
  i=1
  
  #Scores
  plots=list()
  plots_new=list()
  options(ggrepel.max.overlaps = Inf)
  
  for(cell_type in names(cvfit_final_list_maxscore) ){
    
    regression_coef=data.frame(name=rownames(coef(cvfit_final_list_maxscore[[cell_type]], s = "lambda.1se"))[-1], 
                               coef=as.matrix(coef(cvfit_final_list_maxscore[[cell_type]], s = "lambda.1se"))[-1] )
    
    # 
    regression_coef=regression_coef[order(regression_coef$coef, decreasing =TRUE),]
    # 
    regression_coef$order=1:nrow(regression_coef)
    # 
    regression_coef$symbol=representatives$symbol[match(regression_coef$name, representatives$ID)]
    # 
    
    regression_coef$new_motif=representatives$new_motif[match(regression_coef$name, representatives$ID)]
    
    
    writeData(wb, "scores", x=cell_type, startCol=(i-1)*5+1, startRow=1)
    addStyle(wb, "scores", style=createStyle(border="left"), cols=(i-1)*5+1, rows=1)
    writeData(wb, "scores", x=cell_type, startCol=(i-1)*5+2, startRow=1)
    writeData(wb, "scores", x=cell_type, startCol=(i-1)*5+3, startRow=1)
    writeData(wb, "scores", x=cell_type, startCol=(i-1)*5+4, startRow=1)
    writeData(wb, "scores", x=cell_type, startCol=(i-1)*5+5, startRow=1)
    addStyle(wb, "scores", style=createStyle(border="right"), cols=(i-1)*5+5, rows=1)
    
    writeData(wb, "scores", x=regression_coef, startCol=(i-1)*5+1, startRow=2, 
              borders="surrounding", headerStyle = createStyle(border=c("left","right")) )
    
    # Define the style for the rows you want to color (e.g., red background)
    style_red <- createStyle(fgFill = "red")  # FF0000 is the hex code for red
    
    # Apply the style to the desired rows (e.g., rows 2 and 4)
    addStyle(wb, "scores", style = style_red, rows = 2+which(regression_coef$new_motif==TRUE), cols = ((i-1)*5+1):((i-1)*5+5), gridExpand = TRUE)
    
    
    
    
    
    # #Plot regression coefficient as the function of the feature importance
    
    # 
    plots[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point(aes(color = new_motif))+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      #geom_text_repel(data=subset(regression_coef, coef > 0.025, max.overlaps=Inf), aes(label=symbol), box.padding=0.5)+
      
      geom_text_repel(data=subset(regression_coef, (coef > 0.025 | coef < -0.025) & !new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="black", box.padding=0.5, force=10, xlim=c(0, NA), max.time=10)+
      geom_text_repel(data=subset(regression_coef, (coef > 0.025 | coef < -0.025) & new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="red", box.padding=0.5, force=10, xlim=c(0, NA), max.time=10)+
      
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient", color="New motif")+
      theme_minimal()
    
    plots_new[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point(aes(color = new_motif))+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      #geom_text_repel(data=subset(regression_coef, coef > 0.025, max.overlaps=Inf), aes(label=symbol), box.padding=0.5)+
      geom_text_repel(data=subset(regression_coef, (coef > 0.01 ) & new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="red", box.padding=0.5, force=10, xlim=c(0, NA), max.time=10)+
      
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient", color="New motif")+
      theme_minimal()
    
    i=i+1
    
    
    # 
  }
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list_maxscore))))
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs_new_motifs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  # Save the workbook to a file
  saveWorkbook(wb, paste0("/scratch/project_2006203/TFBS/ATAC-seq-peaks/reg_coeff_excels/",cell_type_group,".xlsx"), overwrite = TRUE)
  
#}

