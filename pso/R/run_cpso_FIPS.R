#' Fully Informed Particle Swarm Optimization
#'
#' This function runs the CPSO algorithm to find the minimum of the provided objective function.
#' It initializes a swarm of particles, updates their velocities and positions over iterations,
#' and returns the best solution found.
#'
#' @param objective_function A function to be minimized, which takes a numeric vector as input and returns a numeric value.
#' @param N Integer. The number of particles in the swarm.
#' @param xmin Numeric vector. The lower bounds of the search space.
#' @param xmax Numeric vector. The upper bounds of the search space.
#' @param max_iter Integer. The maximum number of iterations.
#' @param w Numeric. The inertia weight. Defaults to \eqn{\frac{1}{2 \log 2}}.
#' @param c Numeric. The acceleration coefficient. Defaults to \eqn{0.5 + \log 2}.
#' @param rho Numeric. The proportion of the population that acts as informants for each particle. Defaults to \eqn{\varrho = 1 - \Big(1 - \dfrac{1}{N}\Big)^{d}}
#'
#' @return A list containing the best position found and the corresponding fitness value.
#'
#' @details The CPSO algorithm follows the Fully Informed Particle Swarm (FIPS) variant of PSO.
#' The velocity update equation is:
#' \deqn{
#' \mathbf{V}_{t+1}^{(i)} = w \cdot \mathbf{V}_{t}^{(i)} + c \cdot \sum_{j \in \mathcal{I}_i} r_j \cdot (\mathbf{p}_t^{(j)} - \mathbf{X}_t^{(i)})
#' }
#' where:
#' \itemize{
#' \item \eqn{w} is the inertia weight,
#' \item \eqn{c} is the acceleration coefficient,
#' \item \eqn{\mathcal{I}_i} is the set of informants for particle \eqn{i},
#' \item \eqn{\mathbf{p}_t^{(j)}} is the personal best position of the \eqn{j}-th informant,
#' \item \eqn{\mathbf{X}_t^{(i)}} is the current position of the \eqn{i}-th particle.
#' }
#'
#' @seealso \code{\link{initialize_swarm}}, \code{\link{update_velocity}}, \code{\link{update_position}}, \code{\link{evaluate_fitness}}
#'
#' @examples
#' objective_function <- function(x) sum(x^2)
#' result <- run_cpso(objective_function, N = 30, xmin = c(-5, -5), xmax = c(5, 5), max_iter = 100)
#' print(result)
#'
#' @export
run_cpso_FIPS <- function(objective_function,
                          N,
                          xmin,
                          xmax,
                          max_iter,
                          b_max = 20,
                          w = 1 / (2 * log(2)),
                          c = 0.5 + log(2),
                          rho = NULL) {

  # helper functions for discrete kernel
  beta_kernel <- function(N, alpha = 1, beta = 1) {
    x <- seq(0, 1, length.out = N)

    kernel <- dbeta(x, shape1 = alpha, shape2 = beta)

    kernel <- kernel / sum(kernel)

    return(kernel)
  }
  # then a function for building the beta params for iteration t
  # initialize the swarm
  swarm <- initialize_swarm(N, xmin, xmax)
  X <- swarm$X
  V <- swarm$V
  d <- length(xmin)
  # initialize personal bests
  fitness <- evaluate_fitness(X, objective_function)
  personal_bests <- initialize_personal_bests(X, fitness)
  P <- personal_bests$P
  P_fitness <- personal_bests$P_fitness

  # apply default rho value if null
  if(is.null(rho)){
    rho = 1 - (1 - 1/N)^d
  }
  # then declare the number of informants
  m <-  round(rho * N)
  # main FIPS loop
  for (t in 1:max_iter) {
    # get beta for this iteration
    beta_t <- b_max - ((b_max - 1) * (t / max_iter))
    # get kernel for this iteration
    kernel <- beta_kernel(m, alpha = 1, beta = beta_t)
    # sample informants
    informants <- lapply(# for every particle
      1:N, function(i) {
        # remove i from the candidates
        sample(setdiff(1:N, i),
               size = m)  # and sample of size m
      })

    # for every particle...
    # TODO: make this not ass
    for (i in 1:N) {
      # get informant indices for particle i
      I_i <- informants[[i]]
      # get the fitness of the informants
      ordered <- order(P_fitness[I_i])
      inds <- I_i[ordered]
      # informant influence
      I_val <- rowSums(sapply(I_i, function(j) {
        # get the ordered index for this i
        # I_i[ordered[j]]
        # get kernel weight
        phi_j <- kernel[j]
        # random component
        r_j <- runif(d)
        # informant contribution
        r_j * (P[j, ] - X[i, ])
      }))
      # update velocity
      V[i, ] <- w * V[i, ] + c * (P[i, ] - X[i, ]) + c * I_val
    }
    # update positions
    X <- update_position(X, V, xmin, xmax)

    # evaluate fitness
    fitness <- evaluate_fitness(X, objective_function)


    # first, create the probability map from ranks
    probs <- rank_map(fitness)
    # conduct Bernoulli trials with corresponding probability
    inds <- rbinom(N, size = 1, prob = probs) == 1

    # 4.5: Update personal bests
    updated_personal_bests <- update_personal_bests(X, fitness, P, P_fitness)
    P <- updated_personal_bests$P
    P_fitness <- updated_personal_bests$P_fitness

  }

  # determine the global best at the end
  global_best_position <- P[which.min(P_fitness), ]
  global_best_fitness <- min(P_fitness)

  # Step 5: Return the best solution found
  return(list(best_position = global_best_position, best_fitness = global_best_fitness))
}

objective_function <- function(x) sum(x^2)
result <- run_cpso_FIPS(objective_function, N = 30, xmin = c(-5, -5), xmax = c(5, 5), max_iter = 100)
print(result)
