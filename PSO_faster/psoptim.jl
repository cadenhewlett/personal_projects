# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions

function psoptim(par, fn, lower = -1, upper = 1)
    # ------------------- #
    # ------ LOCAL ------ #
    # ------------------- #
    function fn1(par)
        fn(par) / p.fnscale
    end
    return(2)
end    