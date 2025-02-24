#' Chaotic Particle Swarm Optimization with Nonlinear Dynamic Acceleration
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
run_cpso <- function(objective_function,
                     N,
                     xmin,
                     xmax,
                     max_iter,
                     explore_prop = 0.50,
                     w = 1 / (2 * log(2)),
                     c = 0.5 + log(2),
                     rho = NULL,
                     coef_method = "NDAC",
                     mu = 3.99) {
  # initialize the swarm
  swarm <- initialize_swarm(N, xmin, xmax, chaotic = TRUE, mu = mu)
  k = 1 #internal tuning parameter
  X <- swarm$X
  V <- swarm$V
  D <- length(xmin)
  # initialize personal bests
  fitness <- evaluate_fitness(X, objective_function)
  personal_bests <- initialize_personal_bests(X, fitness)
  P <- personal_bests$P
  P_fitness <- personal_bests$P_fitness
  # initial global best - fitness and location
  G <- P[which.min(P_fitness), ]
  G_fitness <- min(P_fitness)
  # initialize acceleration coefficients
  C1 <- initialize_coefs(max_iter, method = coef_method)$c1
  C2 <- initialize_coefs(max_iter, method = coef_method)$c2
  # c3 = c1 * (1 - exp(- c2 * k))
  C3 <- C1 *  (1 - exp( - C2 * k))
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
