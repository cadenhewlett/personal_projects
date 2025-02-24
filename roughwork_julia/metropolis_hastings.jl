# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions
# ------------------------- #
# -- Metropolis Hastings -- #
# ------------------------- #
function metropolis_hastings(gamma::Function, initial_point::Float64, n_iters::Int, proposal_sd = 1)
    samples = Vector{Float64}(undef, n_iters)
    dim = length(initial_point)
    current_point = initial_point
    for i ∈ 1:n_iters
        # propose a new location
        proposal = rand(Normal(current_point, proposal_sd), dim)[dim]
        # compute the MH ratio
        ratio = gamma(proposal) / gamma(current_point) 
        # then compute acceptance Bernoulli
        if rand() < ratio 
            current_point = proposal
        else
            # rejection, nothing to do, i.e. we stay at the current point
            nothing
        end
         # no matter if we reject or not, accumulate one sample
         samples[i] = current_point
    end 
    return samples
end
# --------------------------- #
# -- Test on Beta-Binomial -- #
# --------------------------- #

# prior: beta(α, β)
α = 1
β = 2
# observations are Binomial(n, p)
n  = 3
k  = 3
# build beta-binomial model
function beta_binomial(p::Real)
    if p < 0.0 || p > 1.0
        result = 0.0
    else
        result = pdf(Beta(α, β), p) * pdf(Binomial(n, p), k)
    end
    return result
end

# run MCMC
results = metropolis_hastings( beta_binomial , 0.50, 100000)
# compare posterior mean to theoretical value
println(mean(results))
println((α + k)/(α + β + n))