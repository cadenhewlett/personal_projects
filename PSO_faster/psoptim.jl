# --------------------- #
# -- Adding Packages -- #
# --------------------- #
using Distributions
using LinearAlgebra
using Random

function psoptim(par::Union{Number, Nothing, AbstractVector{<:Number}, Vector{Nothing}},
                 fn::Function;
                 lower::Union{Number, AbstractVector{<:Number}} = -1, 
                 upper::Union{Number, AbstractVector{<:Number}} = 1,
                 kwargs...)
    # ----------------- #
    # -- Basic Setup -- #
    # ----------------- #
    npar = length(par)
    lower = float.([lower[mod1(i, length(lower))] for i in 1:npar])
    upper = float.([upper[mod1(i, length(upper))] for i in 1:npar])
    # ------------------- #
    # --- Local Funcs --- #
    # ------------------- #
    function fn1(par)
        fn(par) / p_fnscale
    end
    # M x N uniform matrix scaled to lower and upper values
    function mnunif(m::Int, n::Int, lower, upper)
        # ensure bounds are column vectors
        lower = isa(lower, Number) ? fill(lower, m) : lower
        upper = isa(upper, Number) ? fill(upper, m) : upper
        # generate m*n matrix scaled to the boundaries
        return rand(m, n) .* (upper .- lower) .+ lower
    end
    # compute Euclidean norm of a numeric vector
    function L2_norm(x::AbstractVector{<:Number})
        return sqrt(sum(x .^ 2))
    end
    # printing utility
    function round_any(x::Any, digits::Int)
        if x === nothing
            out = "nothing"
        elseif x == Inf
            out = "∞"
        elseif x == -Inf
            out = "-∞" 
        else 
            out = round(x, digits=digits)
        end
        return out
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
        :rand_order => true,
        :max_restart => Inf,
        :maxit_stagnate => Inf,
        :trace_stats => false,
        :type => 0
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
    p_trace = controls[:trace] > 0
    # scale factor for function
    p_fnscale = controls[:fnscale] 
    p_type = controls[:type]
    # max evaluations (loop or function)
    p_maxit = Int(controls[:maxit])
    p_maxf = controls[:maxf]
    # absolute and relative tolerances for restarts
    p_abstol = controls[:abstol]
    p_reltol = controls[:reltol]
    # reporting cadence
    p_report = controls[:report]
    # starting swarm size
    p_s = isnothing(controls[:s]) ? Int(floor(10 + sqrt(npar))) : Int(controls[:s])
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
    # compute max velocity
    p_vmax = isnothing(controls[:v_max]) ? controls[:v_max] : controls[:v_max] * p_d
    # random ordering of swarm movement
    p_randorder = controls[:rand_order]
    # maximum restarts and stagnation threshold
    p_maxrestart = controls[:max_restart]
    p_maxstagnate = controls[:maxit_stagnate]
    # collect detailed stats 
    p_trace_stats =  Bool(controls[:trace_stats])
    # scale tolerance by dimension of search space
    if (p_reltol != 0)
        p_reltol = p_reltol*p_d
    end
    # report config before running
    if p_trace
       println("Starting Particle Swarm...")
       @info string("S = ", p_s, ", K = ", controls[:k], ", p = ", round(p_p, digits = 4), 
            ", w0 = ", round(p_w0, digits=4), ", w1 = ", round(p_w1, digits=4), 
            ", c_p = ", round(p_c_p, digits=4), ", c_g = ", round(p_c_g, digits=4),
            )
        @info string(
            "v_max = ", round_any(p_vmax, 4), ", d = ", round_any(p_d, 4),
        )
        @info string(
            "maxit = ", p_maxit, ", maxf = ", round_any(p_maxf, 2), ", abstol = ", 
            round_any(p_abstol, 2), ", maxrestart = ", round_any(p_maxrestart, 2),
            ", maxstagnate = ", round_any(p_maxstagnate, 2)
        )
    end
    # setup performance storage
    if p_trace_stats
        nothing # TODO actually implement iterative storage
    end
  
    # ------------------ #
    # ---- PSO Init ---- #
    # ------------------ #
    # declare the population
    X = mnunif(npar, Int(p_s), lower, upper)
    # replace first column with initial guess if exists
    if !any(isnothing.(par)) && all(par .>= lower) && all(par .<= upper)
        X[:, 1] = par
    end
    # initalization of velocity matrix
    if controls[:type] == 0
        V = mnunif(npar, Int(p_s), lower, upper) - X 
    end
    # if declared, restrict V_0 to maximal velocity
    if !isnothing(p_vmax)
        # temporary object for magnitude of velocity 
        temp = L2_norm.(eachcol(V))
        # rescaling factor wrt maximum
        temp = min.(temp, p_vmax) ./ temp
        # then rescale each particle's V by above factor
        V = V * diagm(temp)
    end
    # initial function evaluation matrix
    f_x = fn1.(eachcol(X))
    # initial function evaluations = swarm size
    stats_feval = p_s
    # initial personal best matrix is X_0
    P = copy(X)
    # similarly, initial best f(x) is f(X_0)
    f_p = copy(f_x)
    # no improvement in the 0th iteration
    P_improved = falses(p_s)
    # extract initial best guess g_0
    i_best = argmin(f_p)
    # get the value of this best guess
    error = f_p[i_best]
    init_links = true
    # initial reporting and storage
    if p_trace && p_report == 1
        println("Iteration 1: fitness = ", round(error, digits = 4))
        if p_trace_stats
            nothing # TODO actually implement iterative storage
        end
    end
    # ----------------------- #
    # ------ Main Loop ------ #
    # ----------------------- #
    stats_iter = 0
    stats_restart = 0
    stats_stagnate = 0
    while  (stats_iter < p_maxit && stats_feval < p_maxf && error > p_abstol && stats_restart < p_maxrestart && stats_stagnate < p_maxstagnate)
        # t <- t + 1
        stats_iter += 1
        # if using informants and links
        if (p_p != 1 && init_links)
            # generate SxS matrix of informant indices
            L = rand(p_s, p_s) .<= p_p
            # ensure that each particle is its own informant
            L[diagind(L)] .= 1
        end
        # generate (un)shuffled indices
        index = p_randorder ? shuffle(1:p_s) : 1:p_s
        # then for each index
        for i ∈ index
            # if using all particles as informants
            if p_p == 1
                # social component is global best
                j = i_best
            else 
                # select the particle index of best informant
                j = findall(L[:, i])[argmin(f_p[L[:, i]])]
            end
            # get t as percentaege completion (wrt maxit)
            t = max(stats_iter / p_maxit, stats_feval / p_maxf)
            # update inertia coefficient if implementing decay
            w_t =  (p_w0 + (p_w1 - p_w0) * t) # TODO: precompute and index with TVAC
            # TODO: implement acceleration updating here  : alternative is to rep p_c_p maxit times and index wrt t
        
            # ------------------------- #
            # ---- Update Velocity ---- #
            # ------------------------- #
            # V_{t + 1}^{(i)} = w_t  V_t^{(i)} + ...
            V[:, i] = w_t .* V[:, i]
            # using SPSO 2007 framework
            if p_type == 0
             # exploration: ... +  r_1c_1 \big( p_t^{(i)} - X^{(i)} \big)
             V[:, i] = V[:, i] + rand(Uniform(0, p_c_p), npar).*(P[:, i] - X[:, i])
             # exploitation: ... +  r_2c_2 \big( g_t - X^{(i)} \big)
             V[:, i] = V[:, i] + rand(Uniform(0, p_c_g), npar).*(P[:, j] - X[:, i])
            end
            # then, check if velocity exceeds cap
            if !isnothing(p_vmax) && abs(p_vmax) != Inf
             # get velocity magnitude
             magV = L2_norm(V[:, i])
             if magV > p_vmax
                # rescale V
                V[:, i] = (p_vmax / temp) .* V[:, i]
             end 
            end
            # ------------------------- #
            # ---- Update Position ---- #
            # ------------------------- #
            # X_{t+1}^{(i)} = X_t^{(i)} + V_t^{(i)}
            X[:, i]  = X[:, i] + V[:, i]
            # check boundaries
            oob = X[:, i] .< lower
            if any(oob)
                # stop particle at the boundaries
                X[oob, i] = lower[oob]
                V[oob, i] .= 0 # bonk!
            end
            # repeat with upper bound
            oob = X[:, i] .> upper
            if any(oob)
                X[oob, i] = upper[oob]
                V[oob, i] .= 0
            end
            # ------------------------ #
            # ---- Update Fitness ---- #
            # ------------------------ #
    
            f_x[i] = fn1(X[:, i])
            stats_feval += 1
            if f_x[i] < f_p[i]
                # update new personal bests & fitness
                P[:, i] = X[:, i]
                # if i == 2  
                #     println("new personal best: Old Best = ", f_p[i], " New Best = ", f_x[i])
                # end
                f_p[i] = f_x[i]
                if f_p[i] < f_p[i_best]
                    # println("Updating i_best: Old Best = ", f_p[i_best], " New Best = ", f_p[i])
                    i_best = i
                end
            end
            # break if we've reached the max number of function evaluations
            if stats_feval >= p_maxf
                break
            end

        end
        # end of swarm movement
        # ------------------ #
        # ---- Restarts ---- #
        #------------------- #
        if p_reltol != 0
            # get distance from best
            d = X .- P[:, i_best]
            # using distance formula, get maximum distance
            d = sqrt(max(sum(d .^ 2, dims=1)))
            # if the largest distance is less than relative tol.
            if d < p_reltol
                # restart 
                X = mnunif(npar, p_s, lower, upper)
                V = (mnunif(npar, p_s, lower, upper) .- X)/2
                # then, check if velocity exceeds cap
                if !isnothing(p_vmax) && abs(p_vmax) != Inf
                    # get velocity magnitude
                    magV = map(L2_norm, eachcol(V))
                    magV = min.(magV, p_vmax) ./ p_vmax
                    # then rescale
                    V = V * diagm(magV)
                end 
                # increase counter 
                stats_restart += 1
                if p_trace 
                    println("It ", stats_iter, ": restarting")
                end
            end
        end
        # ------------------- #
        # ---- Reporting ---- #
        #-------------------- #
        # check if there is improvement
        match_previous = f_p[i_best] == error
        # increment or reset stagnation
        stats_stagnate = match_previous ? stats_stagnate + 1 : 0
        # save fitness in this generation
        error = f_p[i_best]
        # reporting
        if p_trace && (stats_iter % p_report == 0)
            # report stats (diameter if reltol)
            if (p_reltol != 0)
                println("Iteration ", stats_iter, ": fitness = ",
                        round(error, sigdigits = 4), ", swarm diam. =", 
                        round(d, 4)
                    )  
            else
                println("Iteration ", stats_iter, ": fitness = ",
                        round(error, sigdigits = 4))
            end
            # update metrics
            if p_trace_stats
                nothing # TODO actually implement iterative storage
            end
        end
    end
    # end of optimization
    output = (
        par = P[:, i_best], 
        value = f_p[i_best]
    )  
    return (output)
end

f1 = x ->  sum(x.^2)
D = 30
S = 40
# maxF = 2
maxF = 200000
# F1 SETUP
test = psoptim(rand(Uniform(-100, 50), D), f1, 
               lower = fill(-100, D), 
               upper = fill( 100, D), 
               trace = 1, report = Int(maxT/10), maxit = Int(maxF/S), s = S)
 # v_max = 2, 
# print(test)