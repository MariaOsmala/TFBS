#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#or 
#lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

 unload_pkgs <- \(exc=NULL) {
   bpk <- c("splines", "compiler", "tools", "parallel",  "grid", "stats4", "stats", "graphics",  "grDevices", "utils",  "datasets",
"methods",   "base") |> c(exc)
   while (length(setdiff(loadedNamespaces(), bpk)) > 0) {
     lapply(setdiff(loadedNamespaces(), bpk), \(x) {
       try(unloadNamespace(x), silent=FALSE)
     })
   }
 }


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
library("purrr")

#.libPaths("/projappl/project_2006203/project_rpackages_4.2.1/")
.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")
library(ggbreak)

#install.packages("pROC")
#install.packages("PRROC")
library("pROC")
library("PRROC")
cell_line_nro <- as.numeric(commandArgs(trailingOnly = TRUE))

features_type="presence_matrix"
#features_type="max_score_matrix"




#The logistic regression model predicts whether a cCRE i is specific to certain cell type (Class 1) or not (0). 
#The model is trained for all 111 cell types.

#In the simplest model, the covariates are presence/absence (0 or 1)  or counts (0,1,2,...)  
#of motif matches (individual representative PWMS ~1000 features). 
#Or use the binding affinities of motif matches? Positional logistic regression classifier? 
#Use R-packages caret and glmnet. 


setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")
data_path="/scratch/project_2006203/TFBS/"

#representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294 version 1
representatives=read.table("../../PWMs_final_version2.2/metadata_representatives.tsv",sep="\t", header=TRUE) #3933
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

length(unique(rep_motifs))
#1232

#representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")

#motifs=names(representative_motif_matches)
cCREs_list<- readRDS(file = paste0(data_path, "ATAC-seq-peaks/RData/cCREs_list.Rds") ) #
length(cCREs_list)

sum(sapply(cCREs_list, length)) #1 091 805

Adult_Celltypes_mapping <- read_delim("/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#typo
Adult_Celltypes_mapping$`Cell type`[grep("Alverolar Type 2,Immune", Adult_Celltypes_mapping$`Cell type`)]="Alveolar Type 2,Immune"

names(cCREs_list) %in% Adult_Celltypes_mapping$celltype
Adult_Celltypes_mapping$celltype %in% names(cCREs_list) 

#test=data.frame( names(cCREs_list), Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ])

names(cCREs_list)=Adult_Celltypes_mapping$`Cell type`[ match(names(cCREs_list), Adult_Celltypes_mapping$celltype) ]
all_cCREs=unlist(cCREs_list) #which are unique
#length(all_cCREs) #1091805, all are of width 400 bp

#table(width(all_cCREs))

#consider only unique
all_cCREs=unique(all_cCREs)
N=length(all_cCREs) #all possible cCREs #435142



if(features_type=="presence_matrix"){
  presence_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/presence_matrix_version2.2.Rds")) #
  
}else{
  #features_type="max_score_matrix"  
  presence_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/max_score_matrix_version2.2.Rds")) #
}


#Consider only those that are truly unique

cell_type=names(cCREs_list)[cell_line_nro]
#CREs specific to this one

ct_CREs=cCREs_list[[cell_type]]

ol=findOverlaps(ct_CREs, all_cCREs, type="equal", ignore.strand=TRUE) #queryLength: 16481 / subjectLength: 435142

length(unique(ol@from))
length(unique(ol@to))
table(names(all_cCREs[unique(ol@to)]))

#Data frame with rows as samples, columns correspond to the presence of a motif hit, first column as the class label

presence_sparse <- as(presence_matrix, "CsparseMatrix") #435142 x 1232


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


#combined_folds=map2(class1_folds, class2_folds, c)

#other way

results_list=list()
fits_list=list()

random_seed=sample(1:100,1)
print(paste0("seed: ", random_seed))
set.seed(random_seed)

for(ind in 1:ncol(presence_sparse)){
  print(ind)
  class1_folds <- createFolds(which(Class==1), k = 5)
  class2_folds <- createFolds(which(Class==0), k = 5)
  
  class1_folds=lapply(class1_folds,function(x) which(Class==1)[x] )
  class2_folds=lapply(class2_folds,function(x) which(Class==0)[x] )
  
  combined_folds=map2(class1_folds, class2_folds, c)

  
  results <- lapply(combined_folds, function(test_indices) {
    #test_indices=combined_folds[[1]]
    trainData <- presence_sparse[-test_indices, ind]
    testData <- presence_sparse[test_indices, ind]
    
    classes=Class[-test_indices]
    #is a design matrix of dimension n * p, and y is a vector of observations of length n.
    
    data=data.frame(x=as.matrix(trainData), response=classes)
    
    fit=glm(classes ~ x, data=data,  family = binomial()#intercept is included
            )
    
    probabilities=predict(fit, newdata=data.frame(x=testData, response=Class[test_indices]), type="response")
  
    roc_curve <- roc(Class[test_indices], probabilities)
    auc(roc_curve)
    
    pr_curve <- try(pr.curve(scores.class0 = probabilities, weights.class0 = Class[test_indices], curve = TRUE), silent=TRUE)
    
    fit<-list()
    fit[["fit"]]=fit
    fit[["auroc"]]=auc(roc_curve)
    
    if(class(pr_curve)!="try-error"){
      fit[["auprc"]]=pr_curve$auc.integral
    }else{
      fit[["auprc"]]="error"
    }
    
    
    fit
  })
  
  #Fit the whole model
  data=data.frame( x=presence_sparse[, ind], classes=Class)
  fit=glm(classes ~ x, data=data,  family = binomial()#intercept is included
  )

  results_list[[ colnames(presence_matrix)[ind] ]]=results
  fits_list[[ colnames(presence_matrix)[ind] ]]=fit
  
}

if(features_type=="presence_matrix"){
  save.image(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_single_features_version2.2/logistic_regression",cell_type,"_version2.2.RData")) #version 2.2
  
}else{
  #features_type="max_score_matrix" 
  save.image(paste0( data_path, "ATAC-seq-peaks/RData_logistic_regression_single_features_version2.2/logistic_regression_max_score",cell_type,"_version2.2.RData")) #version 2.2
}


