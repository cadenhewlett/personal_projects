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
f2 = x -> ∑(abs.(x)) + ∏(abs.(x)) #  Schwefel's P2.22 function
f3 = function(x) #  Rosenbrock's function
    ∑(d -> (1 - x[d])^2 + 100 * (x[d+1] - x[d]^2)^2, 1:length(x)-1)
end
f4 = function(x) # Noise function 
    U = rand(Uniform(0, 1))
    ∑(d -> (d * x[d]^4), 1:length(x)) + U
end
f5 = function(x) # Schwefel's function
    D = length(x)
    418.9829*D - ∑(d -> (x[d]*sin(sqrt(abs(x[d])))), 1:length(x))
end
f6 = function(x) # Rastrigin's
    ∑(d -> (x[d]^2 - 10*cos(2*π*x[d]) + 10), 1:length(x))
end
f7 = function(x) # Ackley's
    D = length(x)
    T1 = -20*exp(-0.2*sqrt(∑(d -> x[d]^2, 1:length(x)) / D))
    T2 =  -1*exp( ∑(d -> cos(2*π*x[d]), 1:length(x)) / D)
    return T1 + T2 + 20 + exp(1)
end
f8 = function(x) # Griewank's
    1 + ∑(d -> (x[d]^2/4000), 1:length(x)) - ∏(d -> (cos(x[d]/sqrt(d))), 1:length(x))
end
f9 = function(x) # penalized 1
    D = length(x)
    # main penalty term
    u = function(x, a)
       if x > a
        return 100*(x-a)^4
       elseif -a ≤ x ≤ a
        return 0
       else 
        return 100*(-x - a)^4
       end
    end
    # x transform
    y = 1 .+ (x .+ 1)./4
    # main terms
    T1 = 10*(sin(π * y[1])^2)
    T2 = ∑(d -> (y[d] - 1)^2 * (1 + 10*(sin(π * y[d + 1])^2)), 1:(D-1) )
    T3 = (y[D] - 1)^2
    # function value + penalty
    return π*(T1 + T2 + T3)/D + ∑(u.(x, 10))
end
# function list
F = [f1, f2, f3, f4, f5, f6, f7, f8, f9]
# search spaces 
search_L = [-100, -10, -10, -1.28, -500, -5.12, -32, -600, -50]
search_U = [100, 10, 10, 1.28, 500, 5.12, 32, 600, 50]
# initialziation spaces
init_L = [-100, -10, -10, -1.28, -500, -5.12, -32, -600, -50]
init_U = [50, 5, 5, 0.64, 500, 2, 16, 200, 25]

k = 9

result = psoptim(
               rand(Uniform(init_L[k], init_U[k]), D), 
               F[k], 
               lower = fill(search_L[k], D), 
               upper = fill(search_U[k], D), 
               trace = 1, 
               report = Int(maxFE/400), 
               maxit = Int(maxFE/S), 
               s = S)
 # v_max = 2, 
# print(test)