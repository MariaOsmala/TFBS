
library(caret)
library(glmnet)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

data(mtcars)


set.seed(123)
ctrl <- trainControl(method = "cv", number = 10, allowParallel = TRUE)

# Define a grid of hyperparameters. Typically, you might want to test more values.
alpha_vals <- seq(0, 1, by = 0.1)  # Sequence from ridge (0) to lasso (1)
lambda_vals <- 10^seq(-3, 3, length.out = 10)
grid <- expand.grid(alpha = alpha_vals, lambda = lambda_vals)

result <- system.time({ 
model <- train(
  mpg ~ .,
  data = mtcars,
  method = "glmnet",
  trControl = ctrl,
  tuneGrid = grid
)
})

#user  system elapsed 
#0.502   0.019   2.929 

print(model)

stopCluster(cl)


data(BinomialExample)
x <- BinomialExample$x
y <- BinomialExample$y

fit <- glmnet(x, y, family = "binomial", alpha=1) #with lasso, glmnet fits the model for 100 values of lambda by default

#Each curve corresponds to a variable. It shows the path of its coefficient against the l1-norm of the whole coefficient vector as 𝜆
#varies. The axis above indicates the number of nonzero coefficients at the current 𝜆, which is the effective degrees of freedom (df) 
#for the lasso. Users may also wish to annotate the curves: this can be done by setting label = TRUE in the plot command.

plot(fit)

#The number of nonzero coefficients (Df), the percent (of null) deviance explained (%dev) and the value of 𝜆
#(Lambda). Although glmnet fits the model for 100 values of lambda by default, it stops early if %dev does not change sufficiently from one lambda to the next
print(fit)

#We can obtain the model coefficients at one or more 𝜆’s within the range of the sequence:
coef(fit, s = 0.1)

#make predictions at specific 𝜆’s with new input data:
#set.seed(29)
#nx <- matrix(rnorm(5 * 20), 5, 20)
#predict(fit, newx = nx, s = c(0.1, 0.05))

cvfit <- cv.glmnet(x, y,trace.it = TRUE)

#plots the cross-validation curve (red dotted line) along with upper and lower standard deviation curves along the 𝜆
#sequence (error bars). Two special values along the 𝜆
#sequence are indicated by the vertical dotted lines. lambda.min is the value of 𝜆
#that gives minimum mean cross-validated error, while lambda.1se is the value of 𝜆
#that gives the most regularized model such that the cross-validated error is within one standard error of the minimum.
plot(cvfit)


#We can use the following code to get the value of lambda.min and the model coefficients at that value of 𝜆

cvfit$lambda.min
cvfit$lambda.1se

coef(cvfit, s = "lambda.min")

#Note that the coefficients are represented in the sparse matrix format. This is because the solutions along 
#the regularization path are often sparse, and hence it is more efficient in time and space to use a sparse format. 
#If you prefer non-sparse format, pipe the output through as.matrix().

#Predictions can be made based on the fitted cv.glmnet object as well. 
#The code below gives predictions for the new input matrix newx at lambda.min:

predict(cvfit, newx = x[1:5,], s = "lambda.min")