library(Rmpfr)

sine_map <- function(x0, B, r = 0.92, singular = TRUE) {
  x_values <- numeric(B + 1)
  x_values[1] <- x0
  for (i in 1:B) {
    x_values[i + 1] <- r * sin(pi * x_values[i])
  }
  if (singular) {
    return(x_values[B + 1])
  } else {
    return(x_values)
  }
}

dyadic_map <- function(x0, B, singular = TRUE) {
  x_values <- numeric(B + 1)
  x_values[1] <- x0
  for (i in 1:B) {
    x_values[i + 1] <- (2 * x_values[i]) %% 1
  }
  if (singular) {
    return(x_values[B + 1])
  } else {
    return(x_values)
  }
}

chebyshev_map <- function(x0, B, n = 2, singular = TRUE) {
  x_values <- numeric(B + 1)
  x_values[1] <- x0
  for (i in 1:B) {
    x_values[i + 1] <- cos(n * acos(x_values[i]))
  }
  if (singular) {
    return(x_values[B + 1])
  } else {
    return(x_values)
  }
}

gauss_map <- function(x0, B, alpha = 6.20, beta = -0.56, singular = TRUE){
  x_values <- numeric(B + 1)
  x_values[1] <- x0
  for (i in 1:B) {
    x_values[i + 1] <- exp(-alpha* (x_values[i]^2)) + beta
  }
  if (singular) {
    return(x_values[B + 1])
  } else {
    return(x_values)
  }
}


x0 <- abs((sqrt(2) - runif(1))/2)
B = 500
par(mfrow = c(2,3))

results <- matrix(0, nrow = 500, ncol = 3)
for(i in 1:nrow(results)){
  # random x0
  x0 <- runif(1, 0, 1)
  # three diff maps
  results[i, ] <- c(
    sine_map(x0, B = B, singular = FALSE)[B],
    gauss_map(x0, B = B, singular = FALSE)[B],
    chebyshev_map(x0, B = B, singular = FALSE)[B]
  )
}
results
par(mfrow = c(1, 3))
hist(results[, 1], main = "Sine Map x_B Spread")
hist(results[, 2], main = "Gauss Map x_B Spread")
hist(results[, 3], main = "Chebyshev Map x_B Spread")


