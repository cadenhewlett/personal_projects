#' Update Matrix of Personal Bests by Fitness
#'
#' @param X A matrix of size `N x d` representing the current positions of the particles.
#' @param fitness A numeric vector of length `N` representing the current fitness values of the particles.
#' @param P A matrix of size `N x d` representing the current personal best positions of the particles.
#' @param P_fitness A numeric vector of length `N` representing the current personal best fitness values.
#'
#' @return A list containing two elements:
#' \item{P}{An updated matrix of size `N x d` representing the personal best positions of the particles.}
#' \item{P_fitness}{An updated numeric vector of length `N` representing the personal best fitness values.}
#'
#' @seealso \code{\link{initialize_personal_bests}}, \code{\link{evaluate_fitness}}
#'
#' @examples
#' N <- 30
#' swarm <- initialize_swarm(N, c(-5, -5), c(5, 5))
#' fitness <- evaluate_fitness(swarm$X, function(x) sum(x^2))
#' personal_bests <- initialize_personal_bests(swarm$X, fitness)
#' updated_personal_bests <- update_personal_bests(swarm$X, fitness, personal_bests$P, personal_bests$P_fitness)
#' @export
update_personal_bests <- function(X, fitness, P, P_fitness) {
  # input validation
  if (!is.matrix(X) || !is.matrix(P)) {
    stop("X and P must be matrices representing the positions of the particles.")
  }

  if (!is.numeric(fitness) || length(fitness) != nrow(X)) {
    stop("fitness must be a numeric vector with length equal to the number of rows in X.")
  }

  if (!is.numeric(P_fitness) || length(P_fitness) != nrow(P)) {
    stop("`P_fitness` must be a numeric vector with length equal to the number of rows in P.")
  }

  if (nrow(X) != nrow(P) || ncol(X) != ncol(P)) {
    stop("X and P must have the same dimensions.")
  }

  # update personal best positions and fitness values where current fitness is better
  better_mask <- ( fitness < P_fitness )
  P[better_mask, ] <- X[better_mask, ]
  P_fitness[better_mask] <- fitness[better_mask]

  return(list(P = P, P_fitness = P_fitness))
}
