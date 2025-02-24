# example
N <- 150
f <- function(x){return(sum(x^2))}
xmax <- c(5, 5)
xmin <- c(-5, -5)

run_cpso <- function(objective_function,
                     N,
                     xmin,
                     xmax,
                     max_iter,
                     w = 1 / (2 * log(2)),
                     c = 0.5 + log(2),
                     rho = NULL) {
  # initialize the swarm
  swarm <- initialize_swarm(N, xmin, xmax, chaotic = TRUE, mu = 3.90)
  X <- swarm$X
  V <- swarm$V
  D <- length(xmin)
  # initialize personal bests
  fitness <- evaluate_fitness(X, objective_function)
  personal_bests <- initialize_personal_bests(X, fitness)
  P <- personal_bests$P
  P_fitness <- personal_bests$P_fitness
  # initial global best - fitness and location
  G <- X[which.min(P_fitness), ]
  G_fitness <- objective_function(G)
  # initialize acceleration coefficients
  C1 <- initialize_coefs(N, method = "NDAC")$c1
  C2 <- initialize_coefs(N, method = "NDAC")$c2
  # main FIPS loop
  for (t in 1:max_iter) {
    # update positions
    X <- update_position(X, V, xmin, xmax)
    # evaluate fitness
    fitness <- evaluate_fitness(X, objective_function)
    # update personal best
    new_best <- update_personal_bests(X, fitness, P, P_fitness)
    P <- new_best$P
    P_fitness <- new_best$P_fitness
    # update global best
    G <- P[which.min(P_fitness), ]
    G_fitness <- min(P_fitness)
    # update velocity
    for (i in 1:N) {
      # random weights
      r <- runif(D)
      # update velocity
      V[i, ] <- w * V[i, ] + C1[t] * r[1] * (P[i, ] - X[i, ]) + C2[t] * r[2] * (G - X[i, ])
    }
  }

  # determine the global best at the end
  G <- P[which.min(P_fitness), ]
  G_fitness <- min(P_fitness)

  # Step 5: Return the best solution found
  return(list(best_position = G, best_fitness = G_fitness))
}

run_cpso(f, N, max_iter = 30, xmax = xmax, xmin = xmin)$best_position
