# Load required library
library(ggplot2)

# set parameters
mu_min <- 3     # minimum mu value
mu_max <- 4     # maximum mu value
resolution <- 2000
mu_steps <- resolution
iterations <- resolution
last_iterations <- 100  # number of iterations to plot (after transient period)

# Pre-allocate a matrix for efficiency
mu_values <- seq(mu_min, mu_max, length.out = mu_steps)
total_points <- mu_steps * last_iterations
data <- data.frame(mu = numeric(total_points), x = numeric(total_points))

# Counter to track row index
index <- 1
iter <- 1
# iterate over different mu values
for (mu in mu_values) {
  # initial value
  x <- 0.35
  progress <- (iter / mu_steps) * 100
  cat(sprintf("\rProgress = %2.0f%%", progress))
  flush.console()
  # perform iterations
  for (i in 1:iterations) {
    x <- mu * x * (1 - x)
    # store only the last iterations (after transients)
    if (i > (iterations - last_iterations)) {
      data$mu[index] <- mu
      data$x[index] <- x
      index <- index + 1
    }
  }
  iter <- iter + 1
}

# remove unused rows (if any)
data <- data[1:(index-1), ]
# plot using ggplot
p <- ggplot(data, aes(x = mu, y = x)) +
  geom_point(shape = '.', color = 'black') +
  labs(title = "Bifurcation Diagram of Logistic Map",
       x = expression(mu),
       y = expression(x[mu])) +
  geom_vline(xintercept = 3.569945673, color = "red", linetype = "dashed", linewidth = 0.5) +
  theme_bw()
p
# save the plot to a file
ggsave("bifurcation_diagram.png", plot = p, width = 8, height = 6)

