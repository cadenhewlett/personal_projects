using Random
using DataFrames
using Plots
Random.seed!(1928)


# F1 SETUP
S = 5
D = 2
func = function(x)
    return  sum(x.^2)
end

test = psoptim(nothing, 
               func, 
               lower = fill(-100, D), 
               upper = fill( 100, D), 
               trace = 1, 
               report = 1,
               trace_stats = true, 
               maxit = 100,
               s = S)
               
include("psoptim.jl")

plot(log.(test.err))
title!("Iteration vs. Log-Error Plot, SPSO and f_1")
xlabel!("Iteration")
ylabel!("Log(error)")