using LinearAlgebra
using Pkg
using Random
using Distributions
using Base.Threads
using Statistics
using Plots
function mnunif(m::Int, n::Int, lower, upper)
    # ensure bounds are column vectors
    lower = isa(lower, Number) ? fill(lower, m) : lower
    upper = isa(upper, Number) ? fill(upper, m) : upper
    # generate m*n matrix scaled to the boundaries
    return rand(m, n) .* (upper .- lower) .+ lower
end
p_maxit = 10
function w(N, t; γ_max=3.0, γ_min=1e-3, s=u -> 1 - u,)
   u = t / (p_maxit - 1)
   γ = max(γ_max * s(u), γ_min)  # Prevent gamma from hitting 0
   raw = (1:N) .^ (-γ)
   return raw ./ sum(raw)
end

# toy example
p_c_g = 0.5 + log(2)
p_s = 4
p_p = 0.95
npar = 2
X = [1 0 3 -3;
     3 0 3  4]
P = X
f_p = sum(X.^2, dims = 1)
# generate SxS matrix of informant indices
L = rand(p_s, p_s) .<= p_p
# ensure that each particle cannot be its own informant
L[diagind(L)] .= 0
# test i = 1
i = 1
# collect the informants of i
idxs = findall(L[:, i])
# rank them by their performance
R = idxs[sortperm(f_p[idxs]; rev=false)]
# combine weights dimension-wise via (w^{\top} P) 
K_i = P[:, R] * reverse(w(length(R), p_maxit))
U = rand(Uniform(0, p_c_g), npar)
reverse(w(length(R), 1))
(U./2).*(K_i - X[:, i])
(U./2).*(P[:, 1] - X[:, i])

# performance-weighted (descending kernel) for informants
# and also constriction
# survival probability
# if deb, reinsert at a maximal distance from the swarm (based on current diameter)

