using Random
Random.seed!(1928)
include("psoptim.jl")
# ------------------ #
# --- PSO Tests ---- #
# ------------------ #
D = 30
S = 40
maxFE = 200000
reps = 25
# notation 
∑ = sum
∏ = prod
# functions
f1 = x -> ∑(x.^2) # sphere
f2 = x -> ∑(abs.(x)) + ∏(abs.(x)) #  Schwefel's function
f3 = function(x) #  Rosenbrock's function
    ∑(d -> (1 - x[d])^2 + 100 * (x[d+1] - x[d]^2)^2, 1:length(x)-1)
end
f4 = function(x) # Noise function 
    U = rand(Uniform(0, 1), length(x))
    ∑(d -> d * x[d]^4 .+ U[d], 1:length(x))
end
F = [f1, f2, f3]
# search spaces 
search_L = [-100, -10, -10]
search_U = [100, 10, 10]
# initialziation spaces
init_L = [-100, -10, -10]
init_U = [50, 5, 5]

k = 3
init_L[k]
init_U[k]

result = psoptim(
               rand(Uniform(init_L[k], init_U[k]), D), 
               F[k], 
               lower = fill(search_L[k], D), 
               upper = fill(search_U[k], D), 
               trace = 1, 
               report = Int(maxT/10), 
               maxit = Int(maxFE/S), 
               s = S)
 # v_max = 2, 
# print(test)