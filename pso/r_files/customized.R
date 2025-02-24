# objective function
nll_mlr <- function(beta, X, Y, n, k) {
  beta <- matrix(beta, ncol = k, byrow = TRUE)
  logits <- X %*% beta
  exp_logits <- exp(logits)
  probabilities <- exp_logits / rowSums(exp_logits)
  log_prob <- log(probabilities)
  nll <- -sum(Y * log_prob)
  return(nll)
}

# ------------------------------- #
# optimization procedure
optimize_mlr <- function(X, Y, k, initial_beta, iterations = 5) {
  n <- nrow(X)
  p <- ncol(X)
  
  # initialize beta
  beta <- initial_beta
  
  for (i in 1:iterations) {
    # minimize the negative log-likelihood
    fit <- optim(par = beta, fn = nll_mlr, X = X, Y = Y, n = n, k = k, method = "Nelder-Mead")
    # update beta for the next iteration
    beta <- fit$par
  }
  
  # return the final model parameters
  return(matrix(beta, ncol = k, byrow = TRUE))
}
# ------------------------------- #

# example data
X <- model.matrix(~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = train)
Y <- model.matrix(~ Species - 1, data = train)  # Convert factor to one-hot encoding
# initial beta values (for p predictors and k classes)
initial_beta <- rep(0, ncol(X) * 3)  # 3 classes (k = 3)
# optimization
final_beta <- optimize_mlr(X, Y, k = 3, initial_beta = initial_beta, iterations = 5)
# reshape final beta matrix
beta_matrix <- matrix(final_beta, ncol = 3, byrow = TRUE)

# ------------------------------- #

# function to predict class probabilities
predict_probabilities <- function(beta, X, k) {
  logits <- X %*% beta
  exp_logits <- exp(logits)
  probabilities <- exp_logits / rowSums(exp_logits)
  return(probabilities)
}

# predict probabilities for training data
predicted_probs <- predict_probabilities(final_beta, X, k = 3)

# assign predicted labels
predicted_labels <- apply(predicted_probs, 1, which.max)

# compare with actual labels
actual_labels <- apply(Y, 1, which.max)

# compute accuracy
accuracy <- mean(predicted_labels == actual_labels)
cat("Model accuracy: ", round(accuracy, 4), "\n")

# ------------------------------- #