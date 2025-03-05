using LinearAlgebra
using Random
using Distributions
function mnunif(m::Int, n::Int, lower, upper)
    # ensure bounds are column vectors
    lower = isa(lower, Number) ? fill(lower, m) : lower
    upper = isa(upper, Number) ? fill(upper, m) : upper
    # generate m*n matrix scaled to the boundaries
    return rand(m, n) .* (upper .- lower) .+ lower
end
bound = [2, 2]
S = rand(2) .* 4
oob = S .> 2
any( (rand(s, s) .* 4) .> 2 )

S[oob] = bound[oob]
S
