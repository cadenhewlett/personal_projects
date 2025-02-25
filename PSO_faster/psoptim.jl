# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions

function psoptim(par::Union{Number, AbstractVector{<:Number}},
                 fn::Function;
                 lower::Union{Number, AbstractVector{<:Number}} = -1, 
                 upper::Union{Number, AbstractVector{<:Number}} = 1,
                 kwargs...)
    # ----------------- #
    # -- Basic Setup -- #
    # ----------------- #
    npar = length(par)
    lower = float.([lower[mod1(i, length(lower))] for i in 1:npar])
    upper = float.([lower[mod1(i, length(upper))] for i in 1:npar])
    println(lower)
    # ------------------- #
    # --- Local Funcs --- #
    # ------------------- #
    function fn1(par)
        fn(par) / p_fnscale
    end
    function L2_norm(x::AbstractVector{<:Number})
        return sqrt(sum(x .^ 2))
    end

    # ------------------- #
    # -- Default Param -- #
    # ------------------- #
    default_controls = Dict(
        :trace => 0,
        :fnscale => 1,
        :maxit => 1000,
        :maxf => Inf,
        :abstol => -Inf,
        :reltol => 0,
        :report => 10,
        :s => nothing,
        :k => 3,
        :p => nothing,
        :w => 1 / (2 * log(2)),
        :c_p => 0.5 + log(2),
        :c_g => 0.5 + log(2),
        :c_decay => "None",
        :d => nothing,
        :v_max => nothing,
        :rand_order => true
    )
    # extract default control variable names
    con_keys = keys(default_controls)
    # check if we have any illegal keys
    illegal_keys = setdiff(keys(kwargs), con_keys)
    if !isempty(illegal_keys)
        error("Invalid Control Argument(s): ", join(illegal_keys, ", "))
    end
    # check if bounds are of legal size
    if any([lower == -Inf, upper == Inf])
        error("Invalid Upper or Lower Bounds")
    end
    # ------------------ #
    # ---- Controls ---- #
    # ------------------ #
    controls = merge(default_controls, kwargs)
    # extract all control variables
    p_trace = controls[:trace]  
    # scale factor for function
    p_fnscale = controls[:fnscale] 
    # max evaluations (loop or function)
    p_maxit = controls[:maxit]
    p_maxf = controls[:maxf]
    # absolute and relative tolerances for restarts
    p_abstol = controls[:abstol]
    p_reltol = controls[:reltol]
    # reporting cadence
    p_report = controls[:report]
    # starting swarm size
    p_s = isnothing(controls[:s]) ? floor(10 + sqrt(npar)) : controls[:s]
    # average percent of informants
    p_p = isnothing(controls[:p]) ? 1 - (1 - 1 / p_s)^controls[:k] : controls[:p]
    # get the exploitation constant(s)
    p_w = controls[:w]
    # if two constants are supplied, separate them, else duplicate
    p_w0, p_w1 = length(p_w) > 1 ? (p_w[1], p_w[2]) : (p_w, p_w)
    # get the exploration constants
    p_c_p = controls[:c_p] # local 
    p_c_g = controls[:c_g] # global
    # get the region size (default Euclidean norm)
    p_d = isnothing(controls[:d]) ? L2_norm(upper - lower) : controls[:d]

    return (1)
end    
println(psoptim([2], mean, lower=[1,2]))