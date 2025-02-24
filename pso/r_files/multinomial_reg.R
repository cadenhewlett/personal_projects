library(nnet)
library(glmnet)
library(VGAM)
library(mlogit)

# ----- Multinomial Logistic Regression ----- #

set.seed(1928)
# iris (for now)
data(iris)
inds <- sample(1:nrow(iris), size = nrow(iris) * 0.7)
train <- iris[inds, ]
test  <- iris[-inds, ]
# as a factor
train$Species <- as.factor(train$Species)
# train$Petal.Width <- 1 / exp(train$Petal.Width)
# --------------- Neural Net ---------------- #
# fit the multinomial regression model using neural nets
mlr_nnet <- multinom(
  Species ~ Sepal.Length * Sepal.Width * Petal.Length * Petal.Width,
  data = train,
  trace = FALSE,
  reltol = 1.0e-8
)
# ---------------- GLM CV Lasso --------------------#
# prepare the data for glmnet
x <- model.matrix(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
                  data = train)[, -1]
y <- as.factor(train$Species)
# fit MLR model with Lasso (alpha = 1)
mlr_lasso  <- glmnet(x, y, family = "multinomial", alpha = 1)
# cross-validation
cv_lasso <- cv.glmnet(x, y, family = "multinomial", alpha = 1)

# optimal lambda that minimizes the cross-validated error
lambda_lasso <- cv_lasso$lambda.min

# most regularized model such that the error is within one standard error of the minimum
lambda_1se <- cv_lasso$lambda.1se

# ---------------- GLM CV Ridge --------------------#

# fit MLR model with Ridge (alpha = 0)
mlr_ridge  <- glmnet(x, y, family = "multinomial", alpha = 0)
# cross-validation
cv_ridge <- cv.glmnet(x, y, family = "multinomial", alpha = 0)

# optimal lambda that minimizes the cross-validated error
ridge_lambda <- cv_ridge$lambda.min

# ---------------- GLM CV ElasticNet --------------------#

# fit MLR model with Ridge (alpha = 0)
mlr_enet  <- glmnet(x, y, family = "multinomial", alpha = 0.5)
# cross-validation
cv_enet <- cv.glmnet(x, y, family = "multinomial", alpha = 0.5)

# optimal lambda that minimizes the cross-validated error
enet_lmbda <- cv_enet$lambda.min

# ---------------- Comparison --------------------#

pred_nnet <- predict(mlr_nnet, newdata = train, type = "class")
pred_nnet
conf_matrix_nnet <- table(Predicted = pred_nnet, Actual = train$Species)
accuracy_nnet <- sum(diag(conf_matrix_nnet)) / sum(conf_matrix_nnet)

# predictions with lasso from CV lambda
pred_lasso <- predict(mlr_lasso,
                      newx = x,
                      s = lambda_lasso,
                      type = "class")
conf_matrix_lasso <- table(Predicted = pred_lasso, Actual = train$Species)
accuracy_lasso <- sum(diag(conf_matrix_lasso)) / sum(conf_matrix_lasso)

# predictions with ridge from CV lambda
pred_ridge <- predict(mlr_ridge,
                      newx = x,
                      s = ridge_lambda,
                      type = "class")
conf_matrix_ridge <- table(Predicted = pred_ridge, Actual = train$Species)
accuracy_ridge <- sum(diag(conf_matrix_ridge)) / sum(conf_matrix_ridge)

# predictions with enet from CV lambda
pred_enet <- predict(mlr_enet,
                     newx = x,
                     s = enet_lmbda,
                     type = "class")
conf_matrix_enet <- table(Predicted = pred_enet, Actual = train$Species)
accuracy_enet <- sum(diag(conf_matrix_enet)) / sum(conf_matrix_enet)


# ---- Report ---- #
cat("In-Sample Metrics \n")
print(paste("Neural Net", round(accuracy_nnet, 4)))
print(paste("Lasso", round(accuracy_lasso, 4)))
print(paste("Ridge", round(accuracy_ridge, 4)))
print(paste("Elastic Net", round(accuracy_enet, 4)))



# ---------------- Comparison --------------------#
x <- model.matrix(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
                  data = test)[, -1]
pred_nnet <- predict(mlr_nnet, newdata = test, type = "class")

conf_matrix_nnet <- table(Predicted = pred_nnet, Actual = test$Species)
accuracy_nnet <- sum(diag(conf_matrix_nnet)) / sum(conf_matrix_nnet)

# predictions with lasso from CV lambda
pred_lasso <- predict(mlr_lasso,
                      newx = x,
                      s = lambda_lasso,
                      type = "class")
conf_matrix_lasso <- table(Predicted = pred_lasso, Actual = test$Species)
accuracy_lasso <- sum(diag(conf_matrix_lasso)) / sum(conf_matrix_lasso)

# predictions with ridge from CV lambda
pred_ridge <- predict(mlr_ridge,
                      newx = x,
                      s = ridge_lambda,
                      type = "class")
conf_matrix_ridge <- table(Predicted = pred_ridge, Actual = test$Species)
accuracy_ridge <- sum(diag(conf_matrix_ridge)) / sum(conf_matrix_ridge)

# predictions with enet from CV lambda
pred_enet <- predict(mlr_enet,
                     newx = x,
                     s = enet_lmbda,
                     type = "class")
conf_matrix_enet <- table(Predicted = pred_enet, Actual = test$Species)
accuracy_enet <- sum(diag(conf_matrix_enet)) / sum(conf_matrix_enet)

# ---- Report ---- #
cat("Out-of-Sample Metrics \n")
print(paste("Neural Net", round(accuracy_nnet, 4)))
print(paste("Lasso", round(accuracy_lasso, 4)))
print(paste("Ridge", round(accuracy_ridge, 4)))
print(paste("Elastic Net", round(accuracy_enet, 4)))