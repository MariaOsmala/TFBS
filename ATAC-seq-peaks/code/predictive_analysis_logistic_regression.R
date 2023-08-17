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


#The logistic regression model predicts whether a cCRE i is specific to certain cell type (Class 1) or not (0). 
#The model is trained for all 111 cell types.

#In the simplest model, the covariates are presence/absence (0 or 1)  or counts (0,1,2,...)  
#of motif matches (individual representative PWMS ~1000 features). 
#Or use the binding affinities of motif matches? Positional logistic regression classifier? 
#Use R-packages caret and glmnet. 


#Database stuff
setwd("/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/")

#source("../code/FishersExact.R")

data_path="/scratch/project_2006203/TFBS/"
#gtf<-readRDS(paste0(data_path,"RProjects/TFBS/gtf.Rds")) #This is correct genome, not needed here

#Motif metadata with ICs and length
representatives=read.table("../../PWMs_final/metadata_representatives.tsv",sep="\t", header=TRUE) #3294
rep_motifs=representatives$ID[which(representatives$new_representative=="YES")]

length(unique(rep_motifs))
#1031

representative_motif_matches <- readRDS("/scratch/project_2006203/TFBS/ATAC-seq-peaks/RData/representative_motif_matches_human_final.Rds")

motifs=names(representative_motif_matches)
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



presence_matrix=readRDS(file = paste0( data_path, "ATAC-seq-peaks/RData/presence_matrix.Rds")) #
#saveRDS(count_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/count_matrix.Rds")) #
#saveRDS(max_score_matrix,  file = paste0( data_path, "ATAC-seq-peaks/RData/max_score_matrix.Rds")) #

#Consider only those that are truly unique

cell_type=names(cCREs_list)[56]
#CREs specific to this one

ct_CREs=cCREs_list[[cell_type]]

ol=findOverlaps(ct_CREs, all_cCREs, type="equal", ignore.strand=TRUE)

length(unique(ol@from))
length(unique(ol@to))
table(names(all_cCREs[unique(ol@to)]))

#Data frame with rows as samples, columns correspond to the presence of a motif hit, first column as the class label

presence_df=as.data.frame(presence_matrix)
presence_df$Class=0
presence_df$Class[ol@to]=1

presence_df$Class=as.factor(presence_df$Class)

#Binary features, does it make sense to scale them


class1_folds <- createFolds(which(presence_df$Class==1), k = 5)
class2_folds <- createFolds(which(presence_df$Class==0), k = 5)


class1_folds=lapply(class1_folds,function(x) which(presence_df$Class==1)[x] )
class2_folds=lapply(class2_folds,function(x) which(presence_df$Class==0)[x] )

library("purrr")
combined_folds=map2(class1_folds, class2_folds, c)

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

## All subsequent models are then run in parallel
#model <- train(y ~ ., data = training, method = "rf")

set.seed(123)
models <- lapply(folds, function(test_indices) {
  #test_indices=combined_folds[[1]]
  trainData <- presence_df[-test_indices, ]
  testData <- presence_df[-test_indices, ]
  
  lasso_model <- train(
    Class ~ ., data = trainData,
    method = "glmnet",
    family = "binomial",
    trControl = trainControl(method = "cv"),
    tuneGrid = expand.grid(alpha = 1, lambda = seq(1, 3, by = 1)) #alpha=1 means lasso, alpha=0 is ridge regression
  )
  
  best_lambda <- lasso_model$bestTune$lambda
  # Extract coefficients for the best lambda
  coefficients <- coef(lasso_model$finalModel, s = best_lambda)
  print(coefficients)
  
  coefs_matrix <- as.matrix(lasso_model$finalModel$beta)
  
  plot(lasso_model$finalModel, xvar = "lambda", label = TRUE)
  
  lasso_predictions <- predict(lasso_model, newdata = testData)
  confusionMatrix(lasso_predictions, testData$Class)
  # Here you can insert the code to train your model using 'trainData'
  lasso_model
})

## When you are done:
stopCluster(cl)


data(iris)
binary_iris <- iris[iris$Species %in% c('setosa', 'versicolor'), ] #N=100, D=5

binary_iris$Species=droplevels(binary_iris$Species)



colMeans(binary_iris[,1:4])
var(binary_iris[,1:4])

binary_iris[,1:4]=scale(binary_iris[,1:4])


set.seed(123)
trainIndex <- createDataPartition(binary_iris$Species, p = 0.7, list = FALSE)

folds <- createFolds(binary_iris$Species, k = 5)
models <- lapply(folds, function(train_indices) {
  trainData <- binary_iris[train_indices, ]
  # Here you can insert the code to train your model using 'trainData'
  # Return the trained model
})

results <- lapply(1:length(folds), function(i) {
  train_indices <- folds[[i]]
  testData <- binary_iris[-train_indices, ]
})


trainData <- binary_iris[trainIndex, ]
testData  <- binary_iris[-trainIndex, ]

set.seed(123)
lasso_model <- train(
  Species ~ ., data = trainData,
  method = "glmnet",
  family = "binomial",
  trControl = trainControl(method = "cv"),
  tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01, 0.3, by = 0.01)) #alpha means lasso, alpha=0 is ridge regression
)

print(lasso_model)

best_lambda <- lasso_model$bestTune$lambda
# Extract coefficients for the best lambda
coefficients <- coef(lasso_model$finalModel, s = best_lambda)
print(coefficients)

coefs_matrix <- as.matrix(lasso_model$finalModel$beta)

plot(lasso_model$finalModel, xvar = "lambda", label = TRUE)

#Extract Coefficients:
#You can extract the coefficients for a specific value of 
#λ using the coef function. Remember that glmnet computes coefficients for many values of λ by default.

lambda_value <- lasso_model$results$lambda[30] # Choose a lambda value
coefficients <- coef(lasso_model, s = lambda_value)
print(coefficients)


#The alpha parameter controls the elastic-net mixing parameter where alpha = 0 is ridge regression 
#and alpha = 1 is lasso regression. You can try values between 0 and 1 for elastic net.


set.seed(123)
elastic_net_model <- train(
  Species ~ ., data = trainData,
  method = "glmnet",
  trControl = trainControl(method = "cv"),
  family = "binomial",
  tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.01, 0.3, by = 0.01))
)
print(elastic_net_model)

best_lambda <- elastic_net_model$bestTune$lambda
# Extract coefficients for the best lambda
coefficients <- coef(elastic_net_model$finalModel, s = best_lambda)
print(coefficients)

coefs_matrix <- as.matrix(elastic_net_model$finalModel$beta)

plot(elastic_net_model$finalModel, xvar = "lambda", label = TRUE)

#Extract Coefficients:
#You can extract the coefficients for a specific value of 
#λ using the coef function. Remember that glmnet computes coefficients for many values of λ by default.

lambda_value <- lasso_model$results$lambda[30] # Choose a lambda value
coefficients <- coef(lasso_model, s = lambda_value)
print(coefficients)

#Model evaluation

lasso_predictions <- predict(lasso_model, newdata = testData)
elastic_net_predictions <- predict(elastic_net_model, newdata = testData)

confusionMatrix(lasso_predictions, testData$Species)
confusionMatrix(elastic_net_predictions, testData$Species)

#Visualize the Regularization Path:
#This is more of a bonus step, but you can visualize the coefficients' paths across different lambda values using the glmnet package directly.


lasso_glmnet <- glmnet(as.matrix(trainData[,-5]), trainData$Species, alpha = 1)
plot(lasso_glmnet, xvar = "lambda", label = TRUE)


