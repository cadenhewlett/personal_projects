penalized_nll <- function(beta, X, Y, n, k, lambda, penalty = "ridge") {
  beta <- matrix(beta, ncol = k, byrow = TRUE)
  logits <- X %*% beta
  exp_logits <- exp(logits)
  probabilities <- exp_logits / rowSums(exp_logits)
  log_prob <- log(probabilities)
  nll <- -sum(Y * log_prob)
  
  # add the penalty term based on the selected regularization
  if (penalty == "ridge") {
    penalty_term <- lambda * sum(beta^2)  # L2 regularization
  } else if (penalty == "lasso") {
    penalty_term <- lambda * sum(abs(beta))  # L1 regularization
  } else {
    stop("Invalid penalty type")
  }
  
  return(nll + penalty_term)
}

# ---------------- Implementation --------------------- #
# Example data
train
X <- model.matrix(~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = train)
Y <- model.matrix(~ Species - 1, data = train)  # Convert factor to one-hot encoding

# Initial beta values (for p predictors and k classes)
initial_beta <- rep(0, ncol(X) * 3)  # For 3 classes (k = 3)

# Set regularization parameter lambda
lambda <- 0.1

# Minimize the penalized NLL with Ridge regularization
fit_ridge <- optim(par = initial_beta, fn = penalized_nll, X = X, Y = Y, 
                   n = nrow(X), k = 3, lambda = lambda, penalty = "ridge", method = "BFGS")

# Minimize the penalized NLL with Lasso regularization
fit_lasso <- optim(par = initial_beta, fn = penalized_nll, X = X, Y = Y, 
                   n = nrow(X), k = 3, lambda = lambda, penalty = "lasso", method = "BFGS")

# Extract the optimized coefficients
beta_ridge <- matrix(fit_ridge$par, ncol = 3, byrow = TRUE)
beta_lasso <- matrix(fit_lasso$par, ncol = 3, byrow = TRUE)

# print(beta_ridge)
# print(beta_lasso)

# Function to predict class probabilities
predict_probabilities <- function(beta, X, k) {
  logits <- X %*% beta
  exp_logits <- exp(logits)
  probabilities <- exp_logits / rowSums(exp_logits)
  return(probabilities)
}

# Predict probabilities using Ridge regularization coefficients
predicted_probs_ridge <- predict_probabilities(beta_ridge, X, k = 3)
# Assign predicted labels
predicted_labels_ridge <- apply(predicted_probs_ridge, 1, which.max)
# Compute accuracy
actual_labels <- apply(Y, 1, which.max)
accuracy_ridge <- mean(predicted_labels_ridge == actual_labels)

cat("Ridge Regularization Accuracy:", accuracy_ridge, "\n")

# Predict probabilities using Ridge regularization coefficients
predicted_probs_lasso <- predict_probabilities(beta_lasso, X, k = 3)
# Assign predicted labels
predicted_labels_lasso <- apply(predicted_probs_lasso, 1, which.max)
# Compute accuracy
actual_labels <- apply(Y, 1, which.max)
accuracy_lasso <- mean(predicted_labels_lasso == actual_labels)

cat("Ridge Regularization Accuracy:", accuracy_lasso, "\n")
