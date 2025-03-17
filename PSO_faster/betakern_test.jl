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
p_s = 30
p_p = 0.35
 # generate SxS matrix of informant indices
 L = rand(p_s, p_s) .<= p_p
 # ensure that each particle is its own informant
 L[diagind(L)] .= 1


 beta_weights = function(N, B; A = 1)
    K = (x -> pdf(Beta(A, B), x)).(collect(range(0, 1, length=N)))
    return K/sum(K)
 end

 betafaster = function(N, B)
    K = (x -> B*(1-x)^(B-1)).(collect(range(0, 1, length=N)))
    return K/sum(K)
 end
 
vcat(beta_weights(5, 3)', betafaster(5,3)')
# performance-weighted (descending kernel) for informants
# and also constriction
# survival probability
# if deb, reinsert at a maximal distance from the swarm (based on current diameter)

# other ideas...
# add a Bayesian twist? P(approaching optima | trajectory)
#                       P(global optima | other particles)