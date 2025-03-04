using LinearAlgebra
using Random
function mnunif(m::Int, n::Int, lower, upper)
    # ensure bounds are column vectors
    lower = isa(lower, Number) ? fill(lower, m) : lower
    upper = isa(upper, Number) ? fill(upper, m) : upper
    # generate m*n matrix scaled to the boundaries
    return rand(m, n) .* (upper .- lower) .+ lower
end
mnunif(4, 11, [1, 2, 3, 4], [3, 8, 9, 10])
s = 5
p_p = 0.25
links = rand(s, s) .<= p_p
links[diagind(links)] .= 1

index = shuffle(1:s)

