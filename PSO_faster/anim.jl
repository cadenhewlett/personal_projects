using Random
using DataFrames
using Plots
include("psoptim.jl")


# Simple R2, 5S setup 
S = 5
D = 2
# with f1
func = function(x)
    return  sum(x.^2)
end
# fit particle swarm
M = 100
pso_fit = psoptim(nothing, 
               func, 
               lower = fill(-100, D), 
               upper = fill( 100, D), 
               trace = 1, 
               report = 1,
               trace_stats = true, 
               maxit = M,
               s = S)

# ----------------- #
# --- ANIMATION --- #
# ----------------- #

# define the grid for the heatmap.
xs = range(-150, stop=150, length=300)
ys = range(-150, stop=150, length=300)
z = [func((x,y)) for x in xs, y in ys]

# extract swarm 
swarm = pso_fit.X

anim = @animate for i in 1:M
    # plot f as a heatmap
    heatmap(xs, ys, z,
            legend = false,
            xlims = (-150, 150), ylims = (-150, 150),
            colorbar = false,
            title = "Iteration $i")
    # overlay swarm iteration
    scatter!(swarm[1, :, i], swarm[2, :, i],
             legend = false,
             markersize = 3,
             markercolor = :white)
end

gif(anim, "swarm_on_heatmap.gif", fps = 10)