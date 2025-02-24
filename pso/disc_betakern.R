
maxit = 50L
s = p.s = 10
m = 2

library(CPSO)

swarm <- initialize_swarm(p.s, -1, 1)
X <- swarm$X
V <- swarm$V

links <- matrix(0, p.s, p.s)
for (i in 1:p.s) {
  # remove diag. from sample space
  positions <- setdiff(1:p.s, i)
  # random sample of m informants
  informants <- sample(positions, m)

  links[i, informants] <- 1
}

beta_kernel <- function(N, alpha = 1.0, beta = 1) {
  x <- seq(0, 1, length.out = N)

  kernel <- dbeta(x, shape1 = alpha, shape2 = beta)

  kernel <- kernel / sum(kernel)

  return(kernel)
}

b_max = 40

betas <- sapply(1:maxit, function(t) {
  b_max - ((b_max - 1) * (t / maxit))
})


library(ggplot2)
library(reshape2)
df <- melt(links)
colnames(df) <- c("Row", "Col", "Value")

print(ggplot(df, aes(x = Col, y = Row, fill = as.factor(Value))) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "NxN Matrix Visualization", x = "Columns", y = "Rows"))

sapply(betas, function(b){
  beta_kernel(s, beta = b)
})

# plot(1:s, type = "n", ylim = c(0,1), xlab = "x", ylab = "y", main = "Evolution of Discrete Beta Kernels")
#
# for(b in betas){
#   lines(beta_kernel(s, beta  = b), col = rgb(0, 0, 0, alpha = 0.1), type = 'o')
# }
# beta_kernel(s, beta  = 20)
#
# mean(sapply(betas, function(b){
#   beta_kernel(s, beta = b)[1]
# }))
