.libPaths("/projappl/project_2006203/project_rpackages_4.2.1")

library("rtracklayer")
library("dbplyr")
library("dplyr")
#library("motifStack")
#library("universalmotif")
library("ggplot2")
library(ggbreak)
#install.packages("ggbreak")
#library("GRanges")
library("readr")
library("caret")
library("glmnet")
library("reshape2")
library(ggrepel) #
library("gridExtra")

cell_type_group_nro <- as.numeric(commandArgs(trailingOnly = TRUE)) #1:15

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
#1031

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



#for( cell_type_group in names(cell_type_groups) ){ 

  cell_type_group=names(cell_type_groups)[cell_type_group_nro]
  print(cell_type_group)
  aucs=list()
  assessments_list=list()
  
  cvfits_list=list()
  cvfit_final_list=list()
  
  aucs_maxscore=list()
  assessments_list_maxscore=list()
  
  cvfits_list_maxscore=list()
  cvfit_final_list_maxscore=list()
  
  #presence-absence
  for(cell_type in cell_type_groups[[cell_type_group]] ){
    print(cell_type)
    load(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/logistic_regression_",cell_type,"_version2.2.RData"))
    cvfits_list[[cell_type]]=cvfits
    cvfit_final_list[[cell_type]]=cvfit_final
    assessments=list()
  
    for(i in 1:length(cvfits)){
      print(i)
      assessments[[i]]=assess.glmnet(cvfits[[i]], newx = presence_sparse[combined_folds[[i]], ], 
                             newy = Class[combined_folds[[i]]],s = "lambda.1se", family="binomial")
    
    
    }
    assessments_list[[cell_type]]=assessments
    aucs[[cell_type]]=sapply(assessments, function(x) x[["auc"]])
    
  }
  
  #maxscore
  for(cell_type in cell_type_groups[[cell_type_group]] ){
    print(cell_type)
    load(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/logistic_regression_max_score",cell_type,"_version2.2.RData"))
    cvfits_list_maxscore[[cell_type]]=cvfits
    cvfit_final_list_maxscore[[cell_type]]=cvfit_final
    assessments_maxscore=list()
    
    for(i in 1:length(cvfits)){
      print(i)
      assessments_maxscore[[i]]=assess.glmnet(cvfits[[i]], newx = presence_sparse[combined_folds[[i]], ], #the name of the matrix is the same
                                     newy = Class[combined_folds[[i]]],s = "lambda.1se", family="binomial")
      
      
    }
    assessments_list_maxscore[[cell_type]]=assessments_maxscore
    aucs_maxscore[[cell_type]]=sapply(assessments_maxscore, function(x) x[["auc"]])
    
  }
  
  
  # aucs_ctg[[cell_type_group]]=aucs
  # assessments_ctg[[cell_type_group]]=assessments_list
  # cvfits_ctg[[cell_type_group]]=cvfits_list
  # cvfit_final_ctg[[cell_type_group]]=cvfit_final_list
  # 
  # aucs_maxscore_ctg[[cell_type_group]]=aucs_maxscore
  # assessments_maxscore_ctg[[cell_type_group]]=assessments_list_maxscore
  # cvfits_maxscore_ctg[[cell_type_group]]=cvfits_list_maxscore
  # cvfit_maxscore_final_ctg[[cell_type_group]]=cvfit_final_list_maxscore
  
  
  aucs_test=as.data.frame(do.call(cbind, aucs))
  aucs_test$id=1:nrow(aucs_test)
  aucs_test=reshape2::melt(aucs_test,id.vars=c("id"), value.name="AUC", variable.name = "Cell type")
  aucs_test$method="presence"
  
  aucs_test_maxscore=as.data.frame(do.call(cbind, aucs_maxscore))
  aucs_test_maxscore$id=1:nrow(aucs_test_maxscore)
  aucs_test_maxscore=reshape2::melt(aucs_test_maxscore,id.vars=c("id"), value.name="AUC", variable.name = "Cell type")
  aucs_test_maxscore$method="score"
  
  aucs_test=rbind(aucs_test, aucs_test_maxscore)
  
  
  reshape2::melt(aucs_test, )
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_AUC_boxplots/",cell_type_group,".pdf"), width=6*(length(cell_type_groups[[cell_type_group]])/7), height=6 )
  
  gg=ggplot(aucs_test, aes(x=factor(`Cell type`), y=AUC, fill=method)) + 
      geom_boxplot() + 
      labs(title="Cell-type-specific cCRE prediction based on motif matches, cross-validation,\nPerformance of logistic regression on test data",
           x="Cell type ",y="AUC") +
      theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title=cell_type_group)
  plot(gg)
  
  dev.off()

  #Importance of the regression weights
  #Presence
  plots=list()
  plots_new=list()
  plots_top=list()
  plots_top_new=list()
  options(ggrepel.max.overlaps = Inf)
  
  representatives$new_motif=FALSE
  
  representatives$new_motif[representatives$study=="fromYimeng"]=TRUE   
  
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
    
    # Top ten most and least predictive motifs
    
    top_coef <- regression_coef %>%
      arrange(desc(coef)) %>%
      slice(1:10)
    
    bottom_coef <- regression_coef %>%
      arrange(coef) %>%
      slice(1:10)
    
    plots_top[[cell_type]]=ggplot(regression_coef, aes(x=order, y=coef)) +
      geom_point(aes(color=new_motif)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      geom_text_repel(data=top_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      geom_text_repel(data=bottom_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      labs(title=cell_type,
           x="Features sorted by importance", y="Regression coefficient", color="New motif") +
      theme_minimal()
    
    # Top ten most predictive new motifs
    
    top_coef <- regression_coef %>% filter(name %in% (representatives %>% filter(new_motif==TRUE) %>% pull(ID) )) %>%
      arrange(desc(coef)) %>%
      slice(1:10)
    
    plots_top_new[[cell_type]]=ggplot(regression_coef, aes(x=order, y=coef)) +
      geom_point(aes(color=new_motif)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      geom_text_repel(data=top_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      labs(title=cell_type,
           x="Features sorted by importance", y="Regression coefficient", color="New motif") +
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
  }
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list))))
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_new_motifs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_top_motifs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_top, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_new_top_motifs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_top_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  
  
  #Scores
  plots=list()
  plots_new=list()
  plots_top=list()
  plots_top_new=list()
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
    
    
    # Top ten most and least predictive motifs
    
    top_coef <- regression_coef %>%
      arrange(desc(coef)) %>%
      slice(1:10)
    
    bottom_coef <- regression_coef %>%
      arrange(coef) %>%
      slice(1:10)
    
    plots_top[[cell_type]]=ggplot(regression_coef, aes(x=order, y=coef)) +
      geom_point(aes(color=new_motif)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      geom_text_repel(data=top_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      geom_text_repel(data=bottom_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      labs(title=cell_type,
           x="Features sorted by importance", y="Regression coefficient", color="New motif") +
      theme_minimal()
    
    # Top ten most predictive new motifs
    
    top_coef <- regression_coef %>% filter(name %in% (representatives %>% filter(new_motif==TRUE) %>% pull(ID) )) %>%
      arrange(desc(coef)) %>%
      slice(1:10)
    
    plots_top_new[[cell_type]]=ggplot(regression_coef, aes(x=order, y=coef)) +
      geom_point(aes(color=new_motif)) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      geom_text_repel(data=top_coef, aes(label=symbol, color=new_motif), box.padding=0.5) +
      labs(title=cell_type,
           x="Features sorted by importance", y="Regression coefficient", color="New motif") +
      theme_minimal()
    
    
    
    
    
    plots_new[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point(aes(color = new_motif))+
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
      #geom_text_repel(data=subset(regression_coef, coef > 0.025, max.overlaps=Inf), aes(label=symbol), box.padding=0.5)+
      geom_text_repel(data=subset(regression_coef, (coef > 0.01 ) & new_motif, max.overlaps=Inf), 
                      aes(label=symbol),color="red", box.padding=0.5, force=10, xlim=c(0, NA), max.time=10)+
      
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient", color="New motif")+
      theme_minimal()
    
    
    
    
    # 
  }
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list_maxscore))))
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_new_motifs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_top_motifs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_top, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
  pdf(file = paste0(scratch, "Figures_version2.2/cell_group_reg_coeffs_new_top_motifs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots_top_new, ncol = grid_size, top=cell_type_group))  
  dev.off()
  
#}


save.image(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_version2.2/logistic_regression_processed_",cell_type_group,".RData"))


