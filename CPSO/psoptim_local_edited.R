options(digits = 16)
psoptim_local <- function (par,
                           fn,
                           ...,
                           lower = -1,
                           upper = 1,
                           control = list()) {
  # ------------------- #
  # ------ LOCAL ------ #
  # ------------------- #
  fn1 <- function(par) {
    fn(par, ...) / p.fnscale
  }
  # M x N uniform matrix between boundaries
  mrunif <- function(n, m, lower, upper) {
    return(matrix(runif(n * m, 0, 1), nrow = n, ncol = m) * (upper - lower) +
             lower)
  }
  # produce the ell-2 norm of the input
  norm <- function(x) {
    sqrt(sum(x * x))
  }
  # detect number of parameters from input
  npar <- length(par)
  # detect the variable ranges from the input
  lower <- as.double(rep(lower, , npar))
  upper <- as.double(rep(upper, , npar))
  # get the control list default
  con <- list(
    trace = 0,
    fnscale = 1,
    maxit = 1000L,
    maxf = Inf,
    abstol = -Inf,
    reltol = 0,
    REPORT = 10,
    s = NA,
    k = 3,
    p = NA,
    w = 1 / (2 * log(2)),
    c.p = .5 + log(2),
    c.g = .5 + log(2),
    c.decay = "None",
    d = NA,
    v.max = NA,
    rand.order = TRUE,
    max.restart = Inf,
    maxit.stagnate = Inf,
    vectorize = FALSE,
    hybrid = FALSE,
    hybrid.control = NULL,
    trace.stats = FALSE,
    type = "SPSO2007"
  )

  # save the names
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  # check user input
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  # Argument error checks
  if (any(upper == Inf | lower == -Inf))
    stop("fixed bounds must be provided")
  p.type <- pmatch(con[["type"]], c("SPSO2007", "SPSO2011")) - 1
  if (is.na(p.type))
    stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
  # ------------------------- #
  # ------- CONSTANTS ------- #
  # ------------------------- #
  # provide output on progress?
  p.trace <- con[["trace"]] > 0L
  # scale funcion by 1/fnscale
  p.fnscale <- con[["fnscale"]]
  # maximal number of iterations
  p.maxit <- con[["maxit"]]
  # maximal number of function evaluations
  p.maxf <- con[["maxf"]]
  # absolute tolerance for convergence
  p.abstol <- con[["abstol"]]
  # relative minimal tolerance for restarting
  p.reltol <- con[["reltol"]]
  # output every REPORT iterations
  p.report <- as.integer(con[["REPORT"]])
  # swarm size
  p.s <- ifelse(is.na(con[["s"]]),
                ifelse(p.type == 0,
                       floor(10 + 2 * sqrt(npar)), 40),
                con[["s"]])
  # average % of informants
  p.p <- ifelse(is.na(con[["p"]]),
                # (1 - 1 / N)^D
                1 - (1 - 1 / p.s) ^ con[["k"]],
                con[["p"]])
  # exploitation constant
  p.w0 <- con[["w"]]
  # TODO: check -> is this inertia decay
  if (length(p.w0) > 1) {
    p.w1 <- p.w0[2]
    p.w0 <- p.w0[1]
  } else {
    p.w1 <- p.w0
  }
  # local exploitation constant
  p.c.p <- con[["c.p"]]
  # global exploration constant
  p.c.g <- con[["c.g"]]
  # if using non-constant coefs
  update_accel <- con[["c.decay"]] != "None"
  if(con[["c.decay"]] != "None"){
    # use nonlinear dynamic acceleration (NDAC)
    c1i <- 0.15
    c1f <- log(2)  # .5 +
    # time varying personal best acceleration
    c1t = sapply(0:p.maxit, function(t){
      -(c1f - c1i) * (t/ p.maxit)^2 + c1f
    })
    # time varying global best acceleration
    c2t = sapply(0:p.maxit, function(t){
      c1i *(1 - t/ p.maxit)^2 + c1f * (t /  p.maxit)
    })
    # set starting coef values
    p.c.p <- c1t[1]
    p.c.g <- c2t[1]
  }
  # domain diameter
  p.d <- ifelse(is.na(con[["d"]]), norm(upper - lower), con[["d"]])
  # maximal velocity
  p.vmax <- con[["v.max"]] * p.d
  # process particles in random order?
  p.randorder <- as.logical(con[["rand.order"]])
  # maximal number of restarts
  p.maxrestart <- con[["max.restart"]]
  # maximal number of iterations without improvement
  p.maxstagnate <- con[["maxit.stagnate"]]
  # vectorize?
  p.vectorize <- as.logical(con[["vectorize"]])
  # use hybrid model? i.e. quasi-Newtonian
  if (is.character(con[["hybrid"]])) {
    p.hybrid <- pmatch(con[["hybrid"]], c("off", "on", "improved")) - 1
    if (is.na(p.hybrid))
      stop("hybrid should be one of \"off\", \"on\", \"improved\"")
  } else {
    # use local BFGS search by default
    p.hybrid <- as.integer(as.logical(con[["hybrid"]]))
  }
  # control parameters for hybrid optim
  p.hcontrol <- con[["hybrid.control"]]
  if ("fnscale" %in% names(p.hcontrol)){
    # adapt hyrid optim to function scale, if applicable
    p.hcontrol["fnscale"] <- p.hcontrol["fnscale"] * p.fnscale
  }else{
    p.hcontrol["fnscale"] <- p.fnscale
  }
  # collect detailed stats?
  p.trace.stats <- as.logical(con[["trace.stats"]])
  # report before running
  if (p.trace) {
    message(
      "S=", p.s,
    ", K=", con[["k"]],
    ", p=", signif(p.p, 4),
    ", w0=",signif(p.w0, 4),
    ", w1=",signif(p.w1, 4),
    ", c.p=",signif(p.c.p, 4),
    ", c.g=",signif(p.c.g, 4)
    )
    message(
      "v.max=",signif(con[["v.max"]], 4),
      ", d=", signif(p.d, 4),
      ", vectorize=", p.vectorize,
      ", hybrid=", c("off", "on", "improved")[p.hybrid + 1]
    )
    if (p.trace.stats) {
      stats.trace.it <- c()
      stats.trace.error <- c()
      stats.trace.f <- NULL
      stats.trace.x <- NULL
    }
  }
  ## Initialization
  if (p.reltol != 0){
    p.reltol <- p.reltol * p.d
  }
  # if using the vector form
  if (p.vectorize) {
    lowerM <- matrix(lower, nrow = npar, ncol = p.s)
    upperM <- matrix(upper, nrow = npar, ncol = p.s)
  }
  # declare the population
  X <- mrunif(npar, p.s, lower, upper)
  if (!any(is.na(par)) &&
      all(par >= lower) && all(par <= upper)){
    X[, 1] <- par
  }
  # default approach for initializing velocity matrix
  if (p.type == 0) {
    V <- (mrunif(npar, p.s, lower, upper) - X) / 2
  }
  # if no maximal velocity is declared
  if (!is.na(p.vmax)) {
    # scale to maximal velocity
    temp <- apply(V, 2, norm)
    temp <- pmin.int(temp, p.vmax) / temp
    V <- V %*% diag(temp)
  }
  # initial function evaluation
  f.x <- apply(X, 2, fn1)
  # store for reporting
  stats.feval <- p.s
  # set personal bests to X_0
  P <- X
  # set best f(X) to initial eval.
  f.p <- f.x
  # no improvement in zero-th iteration
  P.improved <- rep(FALSE, p.s)
  # g_0 = best of initial values
  i.best <- which.min(f.p)
  # evaluation of g_0
  error <- f.p[i.best]
  init.links <- TRUE
  # report metrics
  if (p.trace && p.report == 1) {
    message("It 1: fitness=", signif(error, 4))
    if (p.trace.stats) {
      stats.trace.it <- c(stats.trace.it, 1)
      stats.trace.error <- c(stats.trace.error, error)
      stats.trace.f <- c(stats.trace.f, list(f.x))
      stats.trace.x <- c(stats.trace.x, list(X))
    }
  }
  # declare local values
  stats.iter <- 1
  stats.restart <- 0
  stats.stagnate <- 0
  # ----------------------- #
  # ------ MAIN LOOP ------ #
  # ----------------------- #
  while (stats.iter < p.maxit &&
         stats.feval < p.maxf && error > p.abstol &&
         stats.restart < p.maxrestart &&
         stats.stagnate < p.maxstagnate) {
    # t <- t + 1
    stats.iter <- stats.iter + 1
    # if there are not 100% informants
    if (p.p != 1 && init.links) {
      # generate boolean matrix of informants that is N x N
      links <- matrix(runif(p.s * p.s, 0, 1) <= p.p, p.s, p.s)
      # and set all diagonal elements of this matrix to 1
      diag(links) <- TRUE
    }
    ## The swarm moves
    if (!p.vectorize) {
      # evaluate swarm elements in random order
      if (p.randorder) {
        # sample the swarm
        index <- sample(p.s)
      }
      # for every particle
      for (i in index) {
        # if 100% informants
        if (p.p == 1){
          # social component is global best
          j <- i.best
        } else{
          # best informant
          j <- which(links[, i])[which.min(f.p[links[, i]])]
        }

        # 
        temp <- (p.w0 + (p.w1 - p.w0) * max(stats.iter / p.maxit, stats.feval /
                                              p.maxf))
        if(update_accel){
          p.c.p <- c1t[stats.iter]
          p.c.g <- c2t[stats.iter]
        }
        # ------------------------- #
        # ---- UPDATE VELOCITY ---- #
        # ------------------------- #
        # V_{t + 1} = w*V_t + ...
        V[, i] <- temp * V[, i]
        # using 2007
        if (p.type == 0) {
          # exploitation (personal component)
          V[, i] <- V[, i] + runif(npar, 0, p.c.p) * (P[, i] - X[, i])
          # if self is not the best particle already
          if (i != j){
            # exploration: ... + r_1c_1 \big( p_t^{(i)} - X^{(i)} \big)
            V[, i] <- V[, i] + runif(npar, 0, p.c.g) * (P[, j] - X[, i])
          }
        }
        # if we have velocity capping
        if (!is.na(p.vmax)) {
          # get velocity magnitude (L2 Norm)
          temp <- norm(V[, i])
          # if this exceeds the cap...
          if (temp > p.vmax){
            # reduce by the % we are over
            V[, i] <- (p.vmax / temp) * V[, i]
          }
        }
        # ------------------------- #
        # ---- UPDATE POSITION ---- #
        # ------------------------- #
        X[, i] <- X[, i] + V[, i]
        ## Check bounds
        temp <- (X[, i] < lower)
        # if any dimension of the particle is outside the bounds
        if (any(temp)) {
          # set to the boundary and *stop* the particle
          X[temp, i] <- lower[temp]
          V[temp, i] <- 0
        }
        # likewise with the upper bound
        temp <- (X[, i] > upper)
        if (any(temp)) {
          X[temp, i] <- upper[temp]
          V[temp, i] <- 0
        }
        # ------------------------ #
        # ---- UPDATE FITNESS ---- #
        # ------------------------ #
        f.x[i] <- fn1(X[, i])
        stats.feval <- stats.feval + 1
        # if the new fitness is better than the personal best
        if (f.x[i] < f.p[i]) {
          # update location and fitness
          P[, i] <- X[, i]
          f.p[i] <- f.x[i]
          # if this particle's fitness is better than the global best
          if (f.p[i] < f.p[i.best]) {
            # update index of global best
            i.best <- i
          }
        }
        # break if we've reached the max number of function evaluations
        if (stats.feval >= p.maxf){
          break
        }
      }
    }
    # ----------------------- #
    # ---- CHECK FITNESS ---- #
    # ----------------------- #

    # if nonzero relative tolerance
    if (p.reltol != 0) {
      # get distance
      d <- X - P[, i.best]
      # using distance formula, get maximum distances
      d <- sqrt(max(colSums(d * d)))
      # if highest distance is less than relative tolerance
      if (d < p.reltol) {
        # restart
        X <- mrunif(npar, p.s, lower, upper)
        # reinitialize and check velocity matrix
        V <- (mrunif(npar, p.s, lower, upper) - X) / 2
        # checking with maximum velocities
        if (!is.na(p.vmax)) {
          temp <- apply(V, 2, norm)
          temp <- pmin.int(temp, p.vmax) / temp
          V <- V %*% diag(temp)
        }
        # increase the restart counter
        stats.restart <- stats.restart + 1
        if (p.trace){
          message("It ", stats.iter, ": restarting")
        }
      }
    }
    # if no overall improvement from last iteration
    init.links <- f.p[i.best] == error
    # iterate the stagnation counter
    stats.stagnate <- ifelse(init.links, stats.stagnate + 1, 0)
    # and save the best fitness for this generation
    error <- f.p[i.best]
    # if we are not breaking the main execution loop
    if (p.trace && stats.iter %% p.report == 0) {
      # report the trace statistics
      if (p.reltol != 0) { # with or without swarm diameter
        message("It ", stats.iter, ": fitness=",
                signif(error, 4),", swarm diam.=", signif(d, 4))
      } else{
        # otherwise, regular reporting. most common approach
        message("It ", stats.iter, ": fitness=", signif(error, 4))
      }
      # and update all metrics
      if (p.trace.stats) {
        stats.trace.it <- c(stats.trace.it, stats.iter)
        stats.trace.error <- c(stats.trace.error, error)
        stats.trace.f <- c(stats.trace.f, list(f.x))
        stats.trace.x <- c(stats.trace.x, list(X))
      }
    }
  }
  # breaking messages
  if (error <= p.abstol) {
    msg <- "Converged"
    msgcode <- 0
  } else if (stats.feval >= p.maxf) {
    msg <- "Maximal number of function evaluations reached"
    msgcode <- 1
  } else if (stats.iter >= p.maxit) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  } else if (stats.restart >= p.maxrestart) {
    msg <- "Maximal number of restarts reached"
    msgcode <- 3
  } else {
    msg <- "Maximal number of iterations without improvement reached"
    msgcode <- 4
  }
  # report braking message
  if (p.trace){ message(msg) }
  # format outputs
  o <- list(
    par = P[, i.best],
    value = f.p[i.best],
    counts = c(
      "function" = stats.feval,
      "iteration" = stats.iter,
      "restarts" = stats.restart
    ),
    convergence = msgcode,
    message = msg
  )
  # add trace stats information if requested
  if (p.trace &&
      p.trace.stats)
    o <- c(o, list(
      stats = list(
        it = stats.trace.it,
        error = stats.trace.error,
        f = stats.trace.f,
        x = stats.trace.x
      )
    ))

  return(o)
}

set.seed(42)
print(paste("Base Fitness:", pso::psoptim(rep(NA, 10), function(x){sum(x^2)},
                   lower = rep(-5, times = 10), upper = rep(5, times = 10),
                   control =list(s = 40, maxit = 1000, 
                                 trace = 0, REPORT = 100))$value))
print(paste("New Fitness:",psoptim_local(rep(NA, 10), function(x){sum(x^2)},
                    lower = rep(-5, times = 10),
                    upper = rep(5, times = 10),
                    control = list(s = 40, maxit = 1000))$value))
D = 30
S = 40
MAXF = 100000
pso::psoptim(runif(D, -100, 50), 
             function(x){sum(x^2)},
             lower = rep(-100, times = D),
             upper = rep( 100, times = D),
             control = list(s = S, maxit = MAXF,
                            REPORT = MAXF/10))


