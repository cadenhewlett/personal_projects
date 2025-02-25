# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions

function psoptim(par, fn::Function, lower = -1, upper = 1; kwargs... ) 
    # ------------------- #
    # -- Default Param -- #
    # ------------------- #
    # NOTE: variables need new names
    default_controls = Dict(
        :trace => 0,
        :fnscale => 1,
        :maxit => 1000,
        :maxf => Inf,
        :abstol => -Inf,
        :reltol => 0,
        :report => 10,
        :s => nothing,
        :k => 3
        :p => nothing,
        :w => 1 / (2 * log(2)),
        :c.p => 0.5 + log(2),
        :c.g => 0.5 + log(2),
        :c.decay => "None",
        :d => nothing,
        :v.max => nothing,
        :rand.order => true,
        :max.restart => Inf,
        :maxit.stagnate => Inf,
        :trace.stats => false,
        :type => "SPSO2007" # for now
    )
    # extract default control variable names
    con_keys = keys(default_controls)
    # check if we have any illegal keys
    illegal_keys = setdiff(keys(kwargs), con_keys)
    if !isempty(illegal_keys)
        error("Invalid Control Argument(s): ", join(illegal_keys, ", "))
    end
    # otherwise merge
    controls = merge(default_controls, kwargs)
    # ------------------- #
    # ------ LOCAL ------ #
    # ------------------- #
    function fn1(par)
        fn(par) / p.fnscale
    end
    return(2)
end    
println(psoptim(2, mean, s = 10))