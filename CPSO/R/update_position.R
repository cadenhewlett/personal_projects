#' Update Particle Positions
#'
#' Updates the position of each particle in the swarm based on its current velocity.
#' The new positions are computed by adding the velocity to the current position.
#' If a position exceeds the boundaries defined by `xmin` and `xmax`, it is clipped.
#'
#' @param X A matrix of size `N x d` representing the current positions of the particles.
#' @param V A matrix of size `N x d` representing the current velocities of the particles.
#' @param xmin Numeric vector of length `d`. The minimum bounds for each dimension.
#' @param xmax Numeric vector of length `d`. The maximum bounds for each dimension.
#'
#' @return A matrix of size `N x d` representing the updated positions of the particles.
#'
#' @seealso \code{\link{initialize_swarm}}, \code{\link{update_velocity}}, \code{\link{evaluate_fitness}}
#'
#' @examples
#' N <- 30
#' d <- 2
#' swarm <- initialize_swarm(N, c(-5, -5), c(5, 5))
#' V <- matrix(runif(N * d, -0.1, 0.1), nrow = N, ncol = d)
#' updated_positions <- update_position(swarm$X, V, c(-5, -5), c(5, 5))
#' @export
update_position <- function(X, V, xmin, xmax) {

  # input validation
  if (!is.matrix(X) || !is.matrix(V)) {
    stop("X and V must be matrices representing the positions and velocities of the particles.")
  }

  if (length(xmin) != ncol(X) || length(xmax) != ncol(X)) {
    stop("`xmin` and `xmax` must have the same length as the number of columns in `X`.")
  }

  if (!is.numeric(xmin) || !is.numeric(xmax)) {
    stop("`xmin` and `xmax` must be numeric vectors.")
  }

  # update positions
  X_new <- X + V

  # clip positions to the bounds
  X_new <- pmax(pmin(X_new, xmax), xmin)

  return(X_new)
}
