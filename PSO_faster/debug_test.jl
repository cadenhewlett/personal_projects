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
T = 1000

p_c_p = 0#0.5 + log(2)
p_c_g = 0#0.5 + log(2)
#  mnunif(npar, Int(p_s), lower, upper)
# ϕ = c_1 + c_2
ϕ = p_c_p + p_c_g
# χ = 2 / |2 - ϕ - √(ϕ^2 - 4ϕ)|
χ = ϕ >= 4 ? 2 / abs(ϕ - 2 + sqrt(ϕ^2 - 4ϕ)) : 1.0


p_c_p = collect(range(0.5, stop=2.5, length=T))
p_c_g = collect(range(0.5, stop=2.5, length=T))

ϕ = p_c_p .+ p_c_g
# ifelse.(ϕ .>= 4, 2 ./ abs.(ϕ .- 2 .- sqrt.(ϕ.^2 .- 4 .* ϕ)), 1.0)
t = collect(range(1, stop=T, length=T))
T = 1000
χ = function (t)
    h = 0.05   # http://boulph.free.fr/Krzysztof/test%20.pdf
    1 - (t / T)^h
end
plot(χ.(t), title="Linear Decay", xlabel="Step", ylabel="Value", lw=1)