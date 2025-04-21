# parameter is known 
theta <- 5
# sample size
n <- 20
# simulation size
M <- 10000
# squared error loss
L <- function(theta, alpha_x){
  return((theta - alpha_x)^2)
}
# DGP is poisson
DGP <- function(theta) {
  return(rpois(n, theta))
}
# estimator
alpha <- function(x){
  return(mean(x))
}
# loss simulation
simulate_risk <- function(theta, M = 100) {
  losses <- numeric(M)
  for (m in 1:M) {
    x <- DGP(theta)
    estimate <- alpha(x)
    losses[m] <- L(theta, estimate)
  }
  mean(losses)
}
# then do one run 
risk_estimate <- simulate_risk(theta, M)

# evaluate minimaxity
# Define grid of theta values in your parameter space
theta_grid <- seq(1, 10, by = 0.5)

# Compute the risk at each theta
risks <- sapply(theta_grid, simulate_risk, M = 10000)

# Max risk among them
max_risk <- max(risks)
print(max_risk)

# Optional: visualize
plot(theta_grid, risks, type = "b", main = "Frequentist Risk over θ", xlab = "θ", ylab = "Risk")
abline(h = max_risk, col = "red", lty = 2)