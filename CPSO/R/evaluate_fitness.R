#' Get Swarm Fitness
#'
#' @param X A matrix of size `N x d` representing the positions of the particles in the swarm.
#' @param fun A function that takes a numeric vector of length `d` as input and returns a numeric value (fitness).
#'
#' @return A numeric vector of length `N` representing the fitness values of the particles.
#'
#' @examples
#' swarm <- initialize_swarm(10, 1, 2)
#' objective_function <- function(x) { sum(x^2) }
#' fitness <- evaluate_fitness(swarm$X, objective_function)
#' @export
evaluate_fitness <- function(X, fun) {
  # input validation
  if (!is.matrix(X)) {
    stop("X must be a matrix representing the positions of the particles.")
  }
  # input must be a function
  if (!is.function(fun)) {
    stop("objective_function must be a valid function.")
  }
  # swarm must be a matrix
  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("X must have non-zero rows and columns.")
  }

  # apply the objective function to each particle's position
  fitness <- apply(X, 1, fun)

  # make sure dimensionality matches
  if (!is.numeric(fitness) || length(fitness) != nrow(X)) {
    stop("The output of `fun` must be a numeric vector of the same length as the number of rows in X.")
  }

  return(fitness)
}
