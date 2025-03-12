using LinearAlgebra
using Pkg
using Random
using Distributions
using Base.Threads
using Plots
function mnunif(m::Int, n::Int, lower, upper)
    # ensure bounds are column vectors
    lower = isa(lower, Number) ? fill(lower, m) : lower
    upper = isa(upper, Number) ? fill(upper, m) : upper
    # generate m*n matrix scaled to the boundaries
    return rand(m, n) .* (upper .- lower) .+ lower
end
#  mnunif(npar, Int(p_s), lower, upper)
p = 3
N = 100
mnunif(p, N, 0, 1)
Random.seed!(1921)
# simple plot
X = mnunif(2, N, 0, 1)
scatter(X[1, :], X[2, : ], title="Random Uniform Scatter", xlabel="x", ylabel="y")
# TODO: Make this more efficient
pairwise_distances = [norm(X[:, i] - X[:, j]) for i in 1:N, j in 1:N]
histogram(mean.(eachcol(pairwise_distances)), 
          xlabel="Average Pairwise Distances", 
          ylabel="Frequency", 
          title="Mean Pairwise Distances for Uniform Sample",
          legend=false,
          color=:skyblue)