#' Initialize Acceleration Coefficients
#'
#' This function generates time-varying acceleration coefficients (c1 and c2)
#' for a given optimization method over a specified number of iterations.
#'
#' @param max_iter Numeric. The maximum number of iterations for the algorithm.
#' @param method Character. The method used for coefficient initialization.
#'   One of "TVAC", "SCAC", "NDAC", or "SBAC". Default is "NDAC".
#'
#' @return A list with two numeric vectors: `c1` (personal best acceleration)
#'   and `c2` (global best acceleration) for each iteration from 0 to `max_iter`.
#'
#' @examples
#' # initialize coefficients using the default method (NDAC)
#' initialize_coefs(100)
#'
#' # initialize coefficients using the TVAC method
#' initialize_coefs(100, method = "TVAC")
#'
#' @export
initialize_coefs <- function(max_iter, method = "NDAC"){
  # input validation on iter counts
  if (!is.numeric(max_iter) ||
      length(max_iter) != 1 ||
      max_iter <= 0 || floor(max_iter) != max_iter) {
    stop("Error: 'max_iter' must be a positive integer.")
  }
  # define all methods
  TVAC <- function(max_iter, c1i = 0.5, c1f = 2.5, c2i = 2.5, c2f = 0.5){
    # time varying personal best acceleration
    c1t = sapply(0:max_iter, function(t){
      (c1i - c1f) * ((max_iter - t) / max_iter) + c1f
    })
    # time varying global best acceleration
    c2t = sapply(0:max_iter, function(t){
      (c2i - c2f) * ((max_iter - t) / max_iter) + c2f
    })
    # return both as a list to be used in the algo
    return(list(c1 = c1t, c2 = c2t))
  }
  SCAC <- function(max_iter, rho = 2, delta = 0.5){
    # time varying personal best acceleration
    c1t = sapply(0:max_iter, function(t){
      rho * sin(((max_iter - t) / max_iter) * (pi/2)) + delta
    })
    # time varying global best acceleration
    c2t = sapply(0:max_iter, function(t){
      rho * cos(((max_iter - t) / max_iter) * (pi/2)) + delta
    })
    # return both as a list to be used in the algo
    return(list(c1 = c1t, c2 = c2t))
  }
  NDAC <- function(max_iter, c1i = 0.5, c1f = 2.5){
    # time varying personal best acceleration
    c1t = sapply(0:max_iter, function(t){
      -(c1f - c1i) * (t/max_iter)^2 + c1f
    })
    # time varying global best acceleration
    c2t = sapply(0:max_iter, function(t){
      c1i *(1 - t/max_iter)^2 + c1f * (t / max_iter)
    })
    # return both as a list to be used in the algo
    return(list(c1 = c1t, c2 = c2t))
  }
  SBAC <- function(max_iter, c1i = 0.5, c1f = 2.5, lambda = 0.0001){
    # time varying personal best acceleration
    c1t = sapply(0:max_iter, function(t){
      (1 / (1 + exp( -lambda *t/max_iter ) ) ) + 1*(c1f - c1i)*((t/max_iter) - 1)^2
    })
    # time varying global best acceleration
    c2t = sapply(0:max_iter, function(t){
      (1 / (1 + exp( (- lambda * t)/max_iter ) ) ) + (c1f - c1i)*((t/max_iter))^2
    })
    # return both as a list to be used in the algo
    return(list(c1 = c1t, c2 = c2t))
  }
  # extract chosen one
  method <- match.arg(method, c("TVAC", "SCAC", "NDAC", "SBAC"))
  # return based on selected method
  result <- switch(method,
                   TVAC = TVAC(max_iter),
                   SCAC = SCAC(max_iter),
                   NDAC = NDAC(max_iter),
                   SBAC = SBAC(max_iter),
                   stop("Unknown method"))

  return(result)

}
