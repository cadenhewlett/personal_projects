library(CPSO)
# helper functions for discrete kernel
beta_kernel <- function(N, alpha = 1, beta = 1) {
  x <- seq(0, 1, length.out = N)

  kernel <- dbeta(x, shape1 = alpha, shape2 = beta)

  kernel <- kernel / sum(kernel)

  return(kernel)
}
N <- 100
obj <- function(x){
  return(sum(x^2))
}
swarm <- initialize_swarm(N, c(-2,0, -2), c(1, 1, 1))

X <- swarm$X
V <- swarm$V
D <- ncol(X)
fitness <- evaluate_fitness(X, objective_function)
personal_bests <- initialize_personal_bests(X, fitness)

P <- personal_bests$P
P_fitness <- personal_bests$P_fitness

rho <- 1 - (1 - 1/N)^D

m <- round(N*rho)

t = 1; max_iter = 300; b_max = 20
# for every particle
informants <- lapply(
  1:N, function(i) {
    # remove i from the candidates
    sample(setdiff(1:N, i),
           size = m)  # and sample of size m
  })

# compute the personal component
personal <- c1*(P[i, ] - X[i, ])
# compute the global component
gt <- P[which.min(P_fitness), ]
global <- c2*(gt - X[i, ])

# get beta for this iteration
beta_t <- b_max - ((b_max - 1) * (t / max_iter))
# get kernel for this iteration
kernel <- beta_kernel(m, alpha = 1, beta = beta_t)
i = 3; c1 = c2 = c3 = 1
I_i <- informants[[i]]
# sort indices by their personal best fitness
I_prime <- I_i[order(P_fitness[I_i])]
# compute the social component
social <- c3*rowSums(sapply(1:length(I_prime), function(ind){
 # get the associated index of the informant
 j <- I_prime[ind]
 # extract the discretized beta kernel weight
 phi_j <- kernel[ind]
 # compute the social component for this informant
 phi_j*(P[j, ] -  X[i, ])
}))
# new velocity
(social + global)/2 + personal
