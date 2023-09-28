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

representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

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
    load(paste0( data_path, "ATAC-seq-peaks/RData/logistic_regression_",cell_type,".RData"))
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
    load(paste0( data_path, "ATAC-seq-peaks/RData/logistic_regression_max_score",cell_type,".RData"))
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
  
  pdf(file = paste0(scratch, "Figures/cell_group_AUC_boxplots/",cell_type_group,".pdf"), width=12*(length(cell_type_groups[[cell_type_group]])/7), height=4 )
  
  gg=ggplot(aucs_test, aes(x=factor(`Cell type`), y=AUC, fill=method)) + 
      geom_boxplot() + 
      labs(title="Cell-type-specific cCRE prediction based on motif matches, cross-validation,\nPerformance of logistic regression on test data",
           x="Cell type ",y="AUC") +
      theme_minimal()
  plot(gg)
  
  dev.off()

  #Importance of the regression weights
  #Presence
  plots=list()
  options(ggrepel.max.overlaps = Inf)
  
  for(cell_type in names(cvfit_final_list) ){
  
    regression_coef=data.frame(name=rownames(coef(cvfit_final_list[[cell_type]], s = "lambda.1se"))[-1], 
                               coef=as.matrix(coef(cvfit_final_list[[cell_type]], s = "lambda.1se"))[-1] )
    
    # 
    regression_coef=regression_coef[order(regression_coef$coef, decreasing =TRUE),]
    # 
    regression_coef$order=1:nrow(regression_coef)
    # 
    regression_coef$symbol=representatives$symbol[match(regression_coef$name, representatives$ID)]
    # 
    # #Plot regression coefficient as the function of the feature importance
    
    # 
    plots[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point()+geom_text_repel(data=subset(regression_coef, coef > 0.25, max.overlaps=Inf), 
                                                                        aes(label=symbol), box.padding=0.5)+
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient")+
      theme_minimal()
  # 
  }
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list))))
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs/",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size))  
  dev.off()
  
  
  #Scores
  plots=list()
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
    # #Plot regression coefficient as the function of the feature importance
    
    # 
    plots[[cell_type]]=ggplot(regression_coef,aes(x=order, y=coef))+geom_point()+geom_text_repel(data=subset(regression_coef, coef > 0.025, max.overlaps=Inf), 
                                                                                                 aes(label=symbol), box.padding=0.5)+
      labs(title=cell_type,
           x="Features sorted by importance",y="Regression coefficient")+
      theme_minimal()
    # 
  }
  
  #grid size
  grid_size=ceiling(sqrt(length(names(cvfit_final_list_maxscore))))
  
  pdf(file = paste0(scratch, "Figures/cell_group_reg_coeffs/scores_",cell_type_group,".pdf"), width=4*grid_size, height=4*grid_size )
  do.call(grid.arrange, c(plots, ncol = grid_size))  
  dev.off()
  
  
  
  
#}


save.image(paste0( data_path, "ATAC-seq-peaks/RData/logistic_regression_processed_",cell_type_group,".RData"))
stop()


cvfit_final 

presence_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/presence_matrix.Rds")) #
#saveRDS(count_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/count_matrix.Rds")) #
#saveRDS(max_score_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/max_score_matrix.Rds")) #

#Consider only those that are truly unique

cell_type=names(cCREs_list)[cell_line_nro]
#CREs specific to this one

ct_CREs=cCREs_list[[cell_type]]

ol=findOverlaps(ct_CREs, all_cCREs, type="equal", ignore.strand=TRUE) #queryLength: 16481 / subjectLength: 435142

length(unique(ol@from))
length(unique(ol@to))
table(names(all_cCREs[unique(ol@to)]))

#Data frame with rows as samples, columns correspond to the presence of a motif hit, first column as the class label

presence_sparse <- as(presence_matrix, "CsparseMatrix")


#the frequency of the most prevalent value over the second most frequent value (called the “frequency ratio’’), 
#which would be near one for well-behaved predictors and very large for highly-unbalanced data and
#the “percent of unique values’’ is the number of unique values divided by the total number of samples (times 100) 
#that approaches zero as the granularity of the data increases
#nzv <- nearZeroVar(presence_matrix, saveMetrics= TRUE)
#nzv[nzv$nzv,][1:10,]

#filteredData <- presence_matrix[, -nzv]

#Identifying correlated predictors

#highlyCorDescr <- findCorrelation(presence_matrix, cutoff = .75)
#filteredData <- filteredData[,-highlyCorDescr]

#comboInfo <- findLinearCombos(ltfrDesign)
#comboInfo


#presence_df=as.data.frame(presence_matrix)
#presence_df$Class=0
#presence_df$Class[ol@to]=1

#presence_df$Class=as.factor(presence_df$Class)

Class=rep(0, nrow(presence_sparse))
Class[ol@to]=1

#Binary features, does it make sense to scale them?

#class1_folds <- createFolds(which(presence_df$Class==1), k = 5)
#class2_folds <- createFolds(which(presence_df$Class==0), k = 5)


#class1_folds=lapply(class1_folds,function(x) which(presence_df$Class==1)[x] )
#class2_folds=lapply(class2_folds,function(x) which(presence_df$Class==0)[x] )

library("purrr")
#combined_folds=map2(class1_folds, class2_folds, c)

#other way
class1_folds <- createFolds(which(Class==1), k = 5)
class2_folds <- createFolds(which(Class==0), k = 5)

class1_folds=lapply(class1_folds,function(x) which(Class==1)[x] )
class2_folds=lapply(class2_folds,function(x) which(Class==0)[x] )

combined_folds=map2(class1_folds, class2_folds, c)

#library(doParallel)
#cl <- makePSOCKcluster(5)
#registerDoParallel(cl)

## All subsequent models are then run in parallel
#model <- train(y ~ ., data = training, method = "rf")


cvfits <- lapply(combined_folds, function(test_indices) {
  #test_indices=combined_folds[[1]]
  trainData <- presence_sparse[-test_indices, ]
  testData <- presence_sparse[test_indices, ]
  
  classes=Class[-test_indices]
  
  #fit <- glmnet(x=subset(trainData, select = -Class), y=trainData$Class, family = "binomial", alpha=1, trace.it=TRUE)
  #fit <- glmnet(x=trainData, y=classes, family = "binomial", alpha=1, trace.it=TRUE)
  
  #By specifying keep = TRUE in the cv.glmnet call, a matrix of prevalidated predictions are stored in the returned output as the fit.preval component.
  #We can then use this component in the call to assess.glmnet
  cvfit <- cv.glmnet(x=trainData, y=classes,family = "binomial", alpha=1,trace.it = TRUE, type.measure="auc", keep=TRUE)
  
  cvfit
  
  #plots the cross-validation curve (red dotted line) along with upper and lower standard deviation curves along the 𝜆
  #sequence (error bars). Two special values along the 𝜆
  #sequence are indicated by the vertical dotted lines. lambda.min is the value of 𝜆
  #that gives minimum mean cross-validated error, while lambda.1se is the value of 𝜆
  #that gives the most regularized model such that the cross-validated error is within one standard error of the minimum.
  
  #plot(cvfit)
  
  #We can use the following code to get the value of lambda.min and the model coefficients at that value of 𝜆
  
  #cvfit$lambda.min
  #cvfit$lambda.1se
  
  #coef(cvfit, s = "lambda.min")
  
  #Each curve corresponds to a variable. It shows the path of its coefficient against the l1-norm of the whole coefficient vector as 𝜆
  #varies. The axis above indicates the number of nonzero coefficients at the current 𝜆, which is the effective degrees of freedom (df) 
  #for the lasso. Users may also wish to annotate the curves: this can be done by setting label = TRUE in the plot command.
  
  #plot(cvfit$glmnet.fit, label=TRUE)
  
  #The number of nonzero coefficients (Df), the percent (of null) deviance explained (%dev) and the value of 𝜆
  #(Lambda). Although glmnet fits the model for 100 values of lambda by default, it stops early if %dev does not change sufficiently from one lambda to the next
  #print(cvfit$glmnet.fit)
  
  #We can obtain the model coefficients at one or more 𝜆’s within the range of the sequence:
  #coef(fit, s = 0.01)
  
  #predict(cvfit, newx = testData, s = "lambda.min")
  
  #assessment=assess.glmnet(cvfit, newx = testData, newy = Class[test_indices],s = "lambda.min", family="binomial")
  
  #rocs <- roc.glmnet(cvfit$fit.preval, newy = classes, s = "lambda.min")
  
  #roc.glmnet returns a list of cross-validated ROC data, one for each model along the path. 
  #The code below demonstrates how one can plot the output. 
  #The first line identifies the lambda value giving the best area under the curve (AUC). 
  #Then we plot all the ROC curves in grey and the “winner” in red.
  
  #best <- cvfit$index["min",]
  #plot(rocs[[best]], type = "l")
  #invisible(sapply(rocs, lines, col="grey"))
  #lines(rocs[[best]], lwd = 2,col = "red")
  
  #for test data
  #rocs_test <- roc.glmnet(cvfit, newx=testData, newy = Class[test_indices], s = "lambda.min", family="binomial")
  
  #The first argument to confusion.glmnet should be a glmnet or cv.glmnet object (from which predictions can be made), 
  #or a matrix/array of predictions, such as the kept "fit.preval" component in the output of a cv.glmnet call with keep = TRUE
  #cnf <- confusion.glmnet(cvfit, newx = testData, newy = Class[test_indices],s = "lambda.min", family="binomial")
  #print(cnf)
  
  #In the code below, we identify and print the one achieving the smallest classification error.
  #best <- cvfit$index["min",]
  #print(cnf[[best]])
  
  
  
})



#Fit the whole model

cvfit_final <- cv.glmnet(x=presence_sparse, y=Class,family = "binomial", alpha=1,trace.it = TRUE, type.measure="auc", keep=TRUE)

save.image(paste0( data_path, "ATAC-seq-peaks/RData/logistic_regression_",cell_type,".RData"))

# stop()
# 
# plot(cvfit_final)
# 
# #We can use the following code to get the value of lambda.min and the model coefficients at that value of 𝜆
# 
# cvfit_final$lambda.min
# cvfit_final$lambda.1se
# 
# regression_coef=data.frame(name=rownames(coef(cvfit, s = "lambda.min"))[-1], coef=as.matrix(coef(cvfit, s = "lambda.min"))[-1] )
# 
# regression_coef=regression_coef[order(regression_coef$coef, decreasing =TRUE),]
# 
# regression_coef$order=1:nrow(regression_coef)
# 
# regression_coef$symbol=representatives$symbol[match(regression_coef$name, representatives$ID)]
# 
# #Plot regression coefficient as the function of the feature importance
# library(ggrepel)
# 
# ggplot(regression_coef,aes(x=order, y=coef))+geom_point()+geom_text_repel(data=subset(regression_coef, coef > 0.25), 
#                                                                    aes(label=symbol), box.padding=0.5)
# 
# stop()
# 
# lasso_model <- train(
#   Class ~ ., data = trainData,
#   method = "glmnet",
#   family = "binomial",
#   trControl = trainControl(method = "cv"),
#   tuneGrid = expand.grid(alpha = 1, lambda = seq(1, 3, by = 1)) #alpha=1 means lasso, alpha=0 is ridge regression
# )
# 
# best_lambda <- lasso_model$bestTune$lambda
# # Extract coefficients for the best lambda
# coefficients <- coef(lasso_model$finalModel, s = best_lambda)
# print(coefficients)
# 
# coefs_matrix <- as.matrix(lasso_model$finalModel$beta)
# 
# plot(lasso_model$finalModel, xvar = "lambda", label = TRUE)
# 
# lasso_predictions <- predict(lasso_model, newdata = testData)
# confusionMatrix(lasso_predictions, testData$Class)
# # Here you can insert the code to train your model using 'trainData'
# lasso_model
# 
# ## When you are done:
# stopCluster(cl)
# 
# 
# data(iris)
# binary_iris <- iris[iris$Species %in% c('setosa', 'versicolor'), ] #N=100, D=5
# 
# binary_iris$Species=droplevels(binary_iris$Species)
# 
# 
# 
# colMeans(binary_iris[,1:4])
# var(binary_iris[,1:4])
# 
# binary_iris[,1:4]=scale(binary_iris[,1:4])
# 
# 
# set.seed(123)
# trainIndex <- createDataPartition(binary_iris$Species, p = 0.7, list = FALSE)
# 
# folds <- createFolds(binary_iris$Species, k = 5)
# models <- lapply(folds, function(train_indices) {
#   trainData <- binary_iris[train_indices, ]
#   # Here you can insert the code to train your model using 'trainData'
#   # Return the trained model
# })
# 
# results <- lapply(1:length(folds), function(i) {
#   train_indices <- folds[[i]]
#   testData <- binary_iris[-train_indices, ]
# })
# 
# 
# trainData <- binary_iris[trainIndex, ]
# testData  <- binary_iris[-trainIndex, ]
# 
# set.seed(123)
# lasso_model <- train(
#   Species ~ ., data = trainData,
#   method = "glmnet",
#   family = "binomial",
#   trControl = trainControl(method = "cv"),
#   tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01, 0.3, by = 0.01)) #alpha means lasso, alpha=0 is ridge regression
# )
# 
# print(lasso_model)
# 
# best_lambda <- lasso_model$bestTune$lambda
# # Extract coefficients for the best lambda
# coefficients <- coef(lasso_model$finalModel, s = best_lambda)
# print(coefficients)
# 
# coefs_matrix <- as.matrix(lasso_model$finalModel$beta)
# 
# plot(lasso_model$finalModel, xvar = "lambda", label = TRUE)
# 
# #Extract Coefficients:
# #You can extract the coefficients for a specific value of 
# #λ using the coef function. Remember that glmnet computes coefficients for many values of λ by default.
# 
# lambda_value <- lasso_model$results$lambda[30] # Choose a lambda value
# coefficients <- coef(lasso_model, s = lambda_value)
# print(coefficients)
# 
# 
# #The alpha parameter controls the elastic-net mixing parameter where alpha = 0 is ridge regression 
# #and alpha = 1 is lasso regression. You can try values between 0 and 1 for elastic net.
# 
# 
# set.seed(123)
# elastic_net_model <- train(
#   Species ~ ., data = trainData,
#   method = "glmnet",
#   trControl = trainControl(method = "cv"),
#   family = "binomial",
#   tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.01, 0.3, by = 0.01))
# )
# print(elastic_net_model)
# 
# best_lambda <- elastic_net_model$bestTune$lambda
# # Extract coefficients for the best lambda
# coefficients <- coef(elastic_net_model$finalModel, s = best_lambda)
# print(coefficients)
# 
# coefs_matrix <- as.matrix(elastic_net_model$finalModel$beta)
# 
# plot(elastic_net_model$finalModel, xvar = "lambda", label = TRUE)
# 
# #Extract Coefficients:
# #You can extract the coefficients for a specific value of 
# #λ using the coef function. Remember that glmnet computes coefficients for many values of λ by default.
# 
# lambda_value <- lasso_model$results$lambda[30] # Choose a lambda value
# coefficients <- coef(lasso_model, s = lambda_value)
# print(coefficients)
# 
# #Model evaluation
# 
# lasso_predictions <- predict(lasso_model, newdata = testData)
# elastic_net_predictions <- predict(elastic_net_model, newdata = testData)
# 
# confusionMatrix(lasso_predictions, testData$Species)
# confusionMatrix(elastic_net_predictions, testData$Species)
# 
# #Visualize the Regularization Path:
# #This is more of a bonus step, but you can visualize the coefficients' paths across different lambda values using the glmnet package directly.
# 
# 
# lasso_glmnet <- glmnet(as.matrix(trainData[,-5]), trainData$Species, alpha = 1)
# plot(lasso_glmnet, xvar = "lambda", label = TRUE)
# 
# 
