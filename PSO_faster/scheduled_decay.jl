using Plots
# Schedule function s(u), where u = t / (T - 1)
# Choose one of the following for γ(t) = γ_max * s(u)
# Linear:        s(u) = 1 - u
#     Steady, uniform flattening over time
# Quadratic:     s(u) = (1 - u)^2
#     Slow decay early, rapid flattening near the end
# Square root:   s(u) = sqrt(1 - u)
#     Fast flattening early, then gradual
# Exponential:   s(u) = exp(-k * u)        # k > 0
#     Rapid decay early, controlled by k
# Cosine:        s(u) = 0.5 * (1 + cos(π * u))
#     Smooth cosine-shaped annealing

# Sigmoid:       s(u) = 1 - 1 / (1 + exp(-k * (u - 0.5)))   # k > 0
#     S-shaped transition centered around midpoint
# Log ramp:      s(u) = log(1 + a * (1 - u)) / log(1 + a)   # a > 0
#     Logarithmic flattening, tunable via a

function powerlaw_schedule(N, T; γ_max=3.0, s=u->1 - u)
    weights = [
        begin
            u = t / (T - 1)
            γ = γ_max * s(u)
            raw = (1:N) .^ (-γ)
            raw ./ sum(raw)
        end
        for t in 0:T-1
    ]
    return reverse(weights)
end
function plot_weight_schedule(N, T; γ_max=3.0, s=u->1-u)
    weights = powerlaw_schedule(N, T, γ_max=γ_max, s=s)
    colors = cgrad(:viridis, T) 
    plt = plot(xlabel="Index (i)", ylabel="Weight",
               title="Power-law Weights Over Time (γ_max = $γ_max)")
    for t in 1:T
        plot!(1:N, weights[t], color=colors[t], alpha=0.4, label=false)
    end

    return plt
end

# plot
plot_weight_schedule(30, 50, γ_max=4.0)

