library(deSolve)
library(ggplot2)
library(reshape2)
# Define the Lotka-Volterra system
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- alpha * x - beta * x * y
    dy <- delta * x * y - gamma * y
    return(list(c(dx, dy)))
  })
}

# Parameters
parameters <- c(alpha = 1.5, beta = 1, gamma = 1, delta = 0.5)

# Initial state values
state <- c(x = 10, y = 5)

# Time steps
time <- seq(0, 50, by = 0.1)

# Solving the equations
out <- ode(y = state, times = time, func = lotka_volterra, parms = parameters)
out <- as.data.frame(out)
out_long <- reshape2::melt(out, id.vars = "time", measure.vars = c("x", "y"),
                           variable.name = "Population", value.name = "Size")

# Static 3D line plot using scatterplot3d with separate lines for prey and predator
# Plot the dynamics in 3D using scatterplot3d
scatterplot3d(out$time, out$x, out$y, type = "l",
              main = "Plot of Lotka Volterra System",
              xlab = "Time", ylab = "Prey (x)", zlab = "Predator (y)",
              color = "blue", lwd = 2)
