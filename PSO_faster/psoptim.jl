# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions

function psoptim(par, fn::Function, lower = -1, upper = 1; kwargs... ) 
    # ------------------- #
    # -- Default Param -- #
    # ------------------- #
    default_controls = Dict(
        :control1 => 100,
        :control2 => 200,
        :control3 => 300
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
# println(psoptim(2, mean, control2 = 10))