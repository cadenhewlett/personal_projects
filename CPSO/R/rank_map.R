#' Rank Mapping Function
#'
#' This function defines a nonlinear mapping from a set of natural numbers \eqn{A = [1, N] \subset  \mathbb{N}}
#' to the real interval \eqn{[0, 1] \subset \mathbb{R}} using a specified steepness parameter \eqn{p}.
#'
#' The function \eqn{f: A \mapsto [0, 1]} is defined such that:
#' \itemize{
#'   \item \eqn{f(1) = 1}
#'   \item \eqn{f(N) = a}, where \eqn{a \in [0, 1)}
#'   \item \eqn{f(x)} is monotonically decreasing from \eqn{1} to \eqn{N}, meaning \eqn{f(1)} is the largest value and \eqn{f(N)} is the smallest.
#' }
#'
#' If the supplied numeric vector `x` is not a set of natural numbers, the function will automatically rank the elements in `x` and proceed with the largest rank as \eqn{N}.
#'
#' The function is given by the following monotonic transformation of a power function.
#' \deqn{f(x) = a + (1 - a) \cdot \left(1 - \left(\dfrac{x - 1}{N - 1}\right)^p\right)}
#'
#' @param x Numeric vector of values. If the values are not naturals, they will be ranked, and the largest rank will become \eqn{N}.
#' @param N Optional. An integer representing the upper bound of the set \eqn{A = [1, N]}. If not provided, \eqn{N} is set to the largest rank of `x`.
#' @param a A numeric value in the interval \eqn{[0, 1)} that specifies the function's value at \eqn{x = N}.
#' @param p A positive numeric value that controls the steepness of the function. Larger values result in a steeper decrease.
#'
#' @return A numeric vector of the same length as `x`, containing the mapped values \eqn{f(x)}.
#'
#' @examples
#' N = 1000
#' swarm <- initialize_swarm(N, -2, 2)
#' objective_function <- function(x) { sum(x^2) }
#' fitness <- evaluate_fitness(swarm$X, objective_function)
#'
#' trials <- rbinom(N, size = 1, prob = rank_map(fitness)) == 1
#'
#'
#' fit_df <- data.frame(
#'   fitness = fitness,
#'   index  = 1:N
#' )
#' plot(x = fit_df$index, y = fit_df$fitness)
#' points(x = fit_df[!trials, "index"], y = fit_df[!trials, "fitness"], col = 'red')
#' @export

rank_map <- function(x, a = 0.25, p = 2, N = NULL) {
  # input validation
  if (!is.numeric(x)) {
    stop("x must be a numeric vector.")
  }
  if (!is.numeric(a) || a < 0 || a >= 1) {
    stop("a must be a numeric value in the interval [0, 1).")
  }
  if (!is.numeric(p) || p <= 0) {
    stop("p must be a positive numeric value.")
  }

  # if N is not provided, rank the values and set N to the largest rank
  if (is.null(N)) {
    ranks <- rank(x)
    N <- max(ranks)
  } else {
    ranks <- x
    if (any(ranks < 1) || any(ranks > N)) {
      stop("When N is provided, all elements of x must be in the range [1, N].")
    }
  }

  # compute the function value at x
  f_x <- a + (1 - a) * (1 - ((ranks - 1) / (N - 1))^p)
  return(f_x)
}

