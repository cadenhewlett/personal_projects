#' Initialize Matrices of Personal Bests
#'
#' @param X A matrix of size `N x d` representing the positions of the particles.
#' @param fitness A numeric vector of length `N` representing the current fitness values of the particles.
#'
#' @return A list containing two elements:
#' \item{P}{A matrix of size `N x d` representing the personal best positions of the particles.}
#' \item{P_fitness}{A numeric vector of length `N` representing the personal best fitness values.}
#'
#' @seealso \code{\link{initialize_swarm}}, \code{\link{evaluate_fitness}}
#'
#' @examples
#' swarm <- initialize_swarm(30, c(-5, -5), c(5, 5))
#' fitness <- evaluate_fitness(swarm$X, function(x) sum(x^2))
#' personal_bests <- initialize_personal_bests(swarm$X, fitness)
#' @export
initialize_personal_bests <- function(X, fitness) {
  # input validation
  if (!is.matrix(X)) {
    stop("X must be a matrix representing the positions of the particles.")
  }
  # recalling that fitness is an output of `evaluate_fitness`
  if (!is.numeric(fitness) || length(fitness) != nrow(X)) {
    stop("fitness must be a numeric vector with length equal to the number of rows in X.")
  }

  # initialize personal best positions and fitness values
  P <- X
  P_fitness <- fitness

  return(list(P = P, P_fitness = P_fitness))
}
