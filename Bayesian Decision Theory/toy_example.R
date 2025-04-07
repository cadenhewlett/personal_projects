# priors
pi_1 <- function(theta){
  if(theta == 0.5){
    return(0.2)
  } else {
    return(0.8)
  }
}
pi_2 <- function(theta){
  return(1/2)
}
Pi <- list(pi_1, pi_2)
# alphas
a_1 <- function(x){
  return(mean(X))
}
a_2 <- function(x){
  # a(0)=0.25, a(1)=0.75.
  epsilon <- 0.1 # small change in epsilon -> bayes preferred
  # 0.25 -> 0.25+(5*epsilon / 10) for freq
  # 0.25 -> 0.25+(epsilon/10) for bayes
  return(0.25 + (0.50 + epsilon)*x) # key piece for equality
}
A <- list(a_1, a_2)
# x, theta, PDF and loss
X <- c(0,1)
Theta <- c(1/2, 0)
pX <- function(x, theta){
  return(dbinom(x, 1, theta))
}
L <- function(theta, a){
  return(abs(theta - a))
}
# ---------- FREQUENTIST RISK TABLE ------------
frequentist_risk <- matrix(NA, nrow=length(Theta), ncol=length(A))
rownames(frequentist_risk) <- paste("theta =", Theta)
colnames(frequentist_risk) <- c("a1", "a2")
for (i in seq_along(Theta)) {
  for (j in seq_along(A)) {
    theta <- Theta[i]
    risk <- 0
    for (x in X) {
      # expected loss (eq 1)
      risk <- risk + L(theta, A[[j]](x)) * pX(x, theta)
    }
    frequentist_risk[i, j] <- risk
  }
}
print("Frequentist Risk Table:")
print(frequentist_risk)
# ---------- BAYES RISK TABLE ------------
bayes_risk <- matrix(NA, nrow=length(Pi), ncol=length(A))
rownames(bayes_risk) <- c("pi_1", "pi_2 (uniform)")
colnames(bayes_risk) <- c("a1", "a2")
for (k in seq_along(Pi)) {
  prior <- Pi[[k]]
  for (j in seq_along(A)) {
    risk_bayes <- 0
    for (i in seq_along(Theta)) {
      # integrated risk, eq 4
      risk_bayes <- risk_bayes + frequentist_risk[i, j] * prior(Theta[i])
    }
    bayes_risk[k, j] <- risk_bayes
  }
}
print("Bayes Risk Table:")
print(bayes_risk)
# ---------- COMPUTE ------------
# minimax, eq 2
frequentist_minimax <- min(apply(frequentist_risk, 2, max))
print(paste("Frequentist Minimax Risk =", frequentist_minimax))
bayes_maximin <- max(apply(bayes_risk, 1, min))
print(paste("Bayesian Maximin (Bayes Risk under least favourable prior) =", bayes_maximin))


