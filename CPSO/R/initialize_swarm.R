#' Initialize Particle Swarm
#'
#' @param N Integer. The number of particles in the swarm.
#' @param xmin Numeric vector of length `d`. The minimum bound for each dimension.
#' @param xmax Numeric vector of length `d`. The maximum bound for each dimension.
#'
#' @return A list containing two elements:
#' \item{X}{A matrix of size `N x d` representing the initial positions of the particles.}
#' \item{V}{A matrix of size `N x d` representing the initial velocities of the particles.}
#'
#' @examples
#' par(mfrow = c(2, 1))
#' plot( initialize_swarm(50, c(2, 2), c(3, 3), mu = 3.78)$X[, 1],type = 'l', ylab="Chaotic" )
#' plot( initialize_swarm(50, c(2, 2), c(3, 3), chaotic = FALSE)$X[, 1], type = 'l', col = 'red',  ylab="Uniform" )
#' @export
initialize_swarm <- function(N, xmin, xmax, chaotic = TRUE, mu = 3.99) {

  # input validation
  if (length(xmin) != length(xmax)) {
    stop("Length of xmin and xmax must be equal.")
  }
  # both min and max bounds must be numeric
  if (!is.numeric(xmin) || !is.numeric(xmax)) {
    stop("xmin and xmax must be numeric vectors.")
  }
  # they must be correctly ordered
  if (any(xmin >= xmax)) {
    stop("All elements of xmin must be less than corresponding elements of xmax.")
  }
  # N has to contain at least one individual
  if (N <= 0 || !is.numeric(N)) {
    stop("N must be a positive integer.")
  }
  # then, infer dimensionality from domains
  D <- length(xmin)
  # initialize with chaotic schema
  if(chaotic == TRUE){
    # helper function for checking fixed and periodic points
    check_sequence <- function(seqn, epsilon = 1e-3, delta = 1e-4){
      # y_j^{(i)} \gets y_j^{(i)} + \delta, \text{ if }
      # \exists y^{\dagger} \in  {Y_{\text{fixed}}} \text{ s.t. } |y_j^{(i)}
      # - y^{\dagger}| < \varepsilon
      # fixed points we want to look out for
      fixed_points <- c(0, 0.25, 0.5, 0.75, 1,  # general fixed points
                        0.799, 0.451,           # period-2 cycle values (mu ~ 3.0)
                        0.879, 0.215, 0.675,    # period-3 cycle values (mu ~ 3.83)
                        0.728, 0.157 )# other small periodic cycles (mu ~ 3.74)
      # compute all diffs
      diffs <- sapply(fixed_points, function(y_dagger) {
        abs(y_dagger - seqn)
      })
      # if any are too close to a fixed/periodic point
      while (any(diffs < epsilon)) {
        # find which are below the threshold
        bads <- unique(which(diffs < epsilon, arr.ind = TRUE)[, "row"])
        # perturb all such values by adding delta
        seqn[bads] <- seqn[bads] + delta
        # rectify to be within the valid range (delta, 1-delta)
        seqn[bads] <- pmin(pmax(seqn[bads], delta), 1 - delta)
        # recompute differences
        diffs <- sapply(fixed_points, function(y_dagger) {
          abs(y_dagger - seqn)
        })
      }
      return(seqn)
    }
    # simple logistic map
    logistic_map <- function(xn){
      # compute and return x_{n+1}
      x_n1 <- mu * xn * (1 - xn)
      return(x_n1)
    }
    # initialize population
    X <- matrix(0, nrow = N, ncol = D)
    # initialize i = 0
    X[1, ] <- sapply(check_sequence(runif(D, min = 0, max = 1)), logistic_map)
    # then, for all remaining rows
    for (i in 2:N) {
      X[i, ] <- sapply(X[(i - 1), ], logistic_map)
      X[i, ] <- check_sequence(X[i, ])
    }
    # then we just return X to the original scale
    for (j in 1:D) {
      X[, j] <- xmin[j] + X[, j] * (xmax[j] - xmin[j])
    }
  }
  # otherwise initialize positions uniformly within the bounds
  else{
    X <- matrix(runif(N * D, min = xmin, max = xmax), nrow = N, ncol = D)
  }
  # initialize velocities (can start with zero or a small random value)
  V <- matrix(runif(N * D, min = -abs(xmax-xmin), max = abs(xmax-xmin)), nrow = N, ncol = D)
  # return X and V
  return(list(X = X, V = V))
}
