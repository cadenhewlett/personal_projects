using Random
using DataFrames
using Plots
include("psoptim.jl")
include("experiment_sim.jl")
use_log = true
# ----------------- #
# ---- RUN PSO ---- #
# ----------------- #
# Simple R2, 5S setup 
S = 10
D = 2
k = 2
# with f1
func = F[k]
# fit particle swarm
M = 250
pso_fit = psoptim( rand(Uniform(init_L[k], init_U[k]), D), 
                    F[k], 
                    lower = fill(search_L[k], D), 
                    upper = fill(search_U[k], D),  
               trace = 1, 
               report = 1,
               trace_stats = true, 
               maxit = M,
               s = S)

# ----------------- #
# --- ANIMATION --- #
# ----------------- #
# define the grid for the heatmap.
xs = range(search_L[k], stop = search_U[k], length=300)
ys = range(search_L[k], stop = search_U[k], length=300)
z = [func((x,y)) for x in xs, y in ys]
logz = -1*log.(z)
zanim = 1*z
# heatmap parameters
cmap = cgrad(:cividis, 256)
# grid parameters
xtick_values = collect(range(search_L[k], stop = search_U[k], length=5))
ytick_values = collect(range(search_L[k], stop = search_U[k], length=5))
# extract swarm 
swarm = pso_fit.X
zval = use_log ? logz : z
# build animation
anim = @animate for i in 1:M
    # plot f as a heatmap
    contour(xs, ys, zval,
        legend = false,
        xlims = (search_L[k],  search_U[k]), 
        ylims = (search_L[k],  search_U[k]),
        xticks = (xtick_values, string.(xtick_values)),
        yticks = (ytick_values, string.(ytick_values)),
        colorbar = true,
        c = cmap,
        title = "SPSO Performance Function $k in Iteration $i",
        subtitle = "Log-Scale Contours of F$k")
    # overlay swarm iteration
    scatter!(swarm[1, :, i], swarm[2, :, i],
             legend = false,
             markersize = 3,
             markercolor = parse(Colorant, "#C9E0FD"))

end
# save animation
gif(anim, string("anims/swarm_on_heatmap_log", use_log, "f", k, ".gif"), fps = 10)