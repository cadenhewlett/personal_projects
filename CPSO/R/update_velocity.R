#' Update Swarm Velocities by FIPS
#'
#' @param X A matrix of size `N x d` representing the current positions of the particles.
#' @param V A matrix of size `N x d` representing the current velocities of the particles.
#' @param P A matrix of size `N x d` representing the personal best positions of the particles.
#' @param informants A list of length `N`, where each element is a vector of indices representing the informants for each particle.
#' @param w Numeric. The inertia weight, which controls the influence of the current velocity.
#' @param c Numeric. The acceleration coefficient, which controls the influence of the personal best and social components.
#'
#' @return A matrix of size `N x d` representing the updated velocities of the particles.
#'
#' @seealso \code{\link{initialize_swarm}}, \code{\link{initialize_personal_bests}}, \code{\link{update_personal_bests}}, \code{\link{update_position}}
#'
#' @examples
#' N <- 30
#' swarm <- initialize_swarm(N, c(-5, -5), c(5, 5))
#' fitness <- evaluate_fitness(swarm$X, function(x) sum(x^2))
#' personal_bests <- initialize_personal_bests(swarm$X, fitness)
#' informants <- lapply(1:N, function(i) sample(setdiff(1:N, i), size = floor(0.5 * N) - 1))
#' updated_velocities <- update_velocity(swarm$X, swarm$V, personal_bests$P, informants)
#' head(updated_velocities )
#' @export
update_velocity <- function(X, V, P,
                            informants = lapply(1:nrow(X), function(i)
                              sample(setdiff(1:nrow(X), i), size = floor(0.5 * nrow(X)) - 1)),
                            w = 1 / (2 * log(2)),
                            c = 0.5 + log(2)) {
  # input validation
  if (!is.matrix(X) || !is.matrix(V) || !is.matrix(P)) {
    stop(
      "X, V, and P must be matrices representing the positions and velocities of the particles."
    )
  }

  if (!is.list(informants) || length(informants) != nrow(X)) {
    stop("`informants` must be a list of length equal to the number of particles (rows in X).")
  }

  if (!is.numeric(w) || w <= 0) {
    stop("`w` must be a positive numeric value representing the inertia weight.")
  }

  if (!is.numeric(c) || c <= 0) {
    stop("`c` must be a positive numeric value representing the acceleration coefficient.")
  }

  N <- nrow(X)
  d <- ncol(X)

  # Update velocities
  for (i in 1:N) {
    # get informant indices for particle i
    I_i <- informants[[i]]
    # informant influence
    social_component <- rowSums(sapply(I_i, function(j) {
      # random component
      r_j <- runif(d)
      # informant contribution
      r_j * (P[j, ] - X[i, ])
    }))
    # update velocity
    V[i, ] <- w * V[i, ] + c * (P[i, ] - X[i, ]) + c * social_component / length(I_i)
  }

  return(V)
}

