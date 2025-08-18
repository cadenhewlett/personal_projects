library(Rcpp)
library(pso)
library(suncalc)
library(pracma)
library(MASS)
# ============================================================ #
# ||                                                       ||  #
# ||         GENERALIZED KERNEL REGRESSION                 ||  #
# ||     Main Source Code and Helper Functions             ||  #
# ||                                                       ||  #
# ============================================================ #

# ===================== #
# === Global Enums  === #
# ===================== #
r.class <- c(base = 0L, pearson = 1L, deviance = 2L)
g.class <- c(identity =  0L, log = 11, inverse = 2L)

# ============================== #
# === Model Fitting Function === #
# ============================== #
fit_GKR <- function(xt, yt, configs, report = FALSE){
  # ------------------- #
  # -- Configuration -- #
  # ------------------- #
  # configure initial values for optimization
  beta_inits <- apply(xt, MARGIN = 2, function(x) {
    mean(yt, na.rm = T) / (2 * mean(x, na.rm = T))
  })
  # configure ranges
  R <- rbind(
    # kernel parameters
    do.call(rbind, lapply(1:configs$k, function(i) {
      matrix(c(
        # lower
        1e-3, # beta_inits[i] - 10 * abs(beta_inits[i]),
        -3,  -3,
        # upper
        1e2, #beta_inits[i] + 10 * abs(beta_inits[i]),
        3,   3
      ), ncol = 2, byrow = FALSE)
    })),
    # remaining pars
    matrix(c(
      # set up ARMA parameter ranges
      rep(c(-0.5, 0.5), times = configs$p),
      rep(c(-0.5, 0.5), times = configs$q),
      0.01, 0.01
    ),
    ncol = 2,
    byrow = TRUE)
  )
  # infer lower and upper bounds from matrix
  lower <- R[, 1]
  upper <- R[, 2]
  
  # ------------------ #
  # -- Optimization -- #
  # ------------------ #
  
  # optimization target function
  f <- function(par){
    error(xt = xt, 
          yt = yt, 
          params = par,
          configs = configs
    )
  }
  # set up PSO
  pso_fit <- psoptim(
    par = rep(NA, nrow(R)),
    fn  = function(x) {
      -f(x)
    },
    lower = R[, 1],
    upper = R[, 2],
    control = list(
      trace = 100*report, 
      type = "SPSO2007",
      s = 10*configs$k,
      hybrid = FALSE,
      maxit = 400
    )
  )
  # compute fit
  fit <- predict_response_GKR(xt, yt, pso_fit$par, configs)
  # compute shape estimate
  alpha <- compute_shape(yt = yt, mu_t = fit$mu_hat, params = pso_fit$par)
  pso_fit$par[length(pso_fit$par)] <- alpha
  # return parameters and fits
  return(list(
    par = pso_fit$par,
    mu_t = fit$mu_hat,
    yt_tilde = fit$yt_tilde,
    alpha_hat = alpha,
    configs = configs
  ))
}

# =========================== #
# === Prediction Function === #
# =========================== #
predict_response_GKR <- function(xt, yt, params, configs) {
  # default config is no AR, no MA
  configs_default <- list(
    k = max(1, length(params)/3 - 1), 
    p = 0,  q = 0, g = I,
    kernel_type = "gamma",
    use_log = FALSE 
  )
  # modify by user configuration
  configs <- modifyList(configs_default, configs)
  # then extract relevant hyper-parameters
  k <- configs$k
  p <- configs$p
  q <- configs$q
  # ----------- #
  # --- DKR --- #
  # ----------- #
  pars <- as.data.frame(matrix(
    params[1:(3 * k)],
    nrow = k,
    ncol = 3,
    byrow = TRUE
  ))
  colnames(pars) <- c("beta", "logdelta", "logsigma")
  # build kernels for each set of parameters
  kernels <- sapply(1:k, function(i)
    build_kernel(pars$logdelta[i], pars$logsigma[i], type = configs$kernel_type[i]))
  # convert kernel to list if there is only one kernel
  if (k == 1){
    kernels <- list(kernels)}
  # convolve xt with kernels
  xt_conv <- sapply(
    1:k,
    function(i) {
      convolve_kernel(xt[, i], kernels[[i]])
    }, 
    simplify = "matrix"
  )
  # compute kernel contributions
  yt_tilde       <- rowSums(sweep(xt_conv, 2, pars$beta, `*`))
  # ------------------ #
  # --- ARMA Error --- #
  # ------------------ #
  # check for AR params
  phi <- if (p > 0) params[(3*k + 1):(3*k + p)] else 0
  # ibid for MA params
  theta <- if (q > 0) params[(length(params) - q):(length(params) - 1)] else 0
  # then alpha is the last entry (used in log-likelihood)
  # alpha <- compute_shape(yt = yt, mu_t = fit$mu_hat, params = pso_fit$par)
  # use our mu_t computation function
  mu_t <- compute_mu_vec(y = yt, y_tilde = yt_tilde, phi, theta)
  # return results
  list(
    mu_hat         = mu_t,
    yt_tilde       = yt_tilde,
    tau_t          = mu_t - yt_tilde,
    beta           = beta
  )
}
  

# ================================= #
# ======= Density Function  ======= #
# == Of Mean-Parameterized Gamma == #
# ================================= #
dmeangamma <- function(yt, shape, mu = shape, log = FALSE){
  # input validation
  if(any(is.na(mu) | mu < 0) ){
    if(log == FALSE){
      return(0)
    }else{
      return(-1e10)
    }
  }
  if (any(shape < 0)) {
    stop("Shape Parameter must be Positive!")
  } 
  if (any(shape == 0)){
    if(log == FALSE){
      return(0)
    }else{
      return(-1e10)
    }
  }
  if (any((mu / shape) <= 0)) {
    if(log == FALSE){
      return(0)
    }else{
      return(-1e10)
    }
  } 
  # compute log density
  log_density <- ( shape*(log(shape) + log(yt) - log(mu)) - log(yt)
                   - lgamma(shape) - (shape*yt)/mu )
  
  # if yt not in supp(Yt)
  log_density[yt < 0] <- -1e10
  
  # return
  if(log == FALSE){
    return(exp(log_density))
  } else{
    return(log_density)
  }
}

# ================================= #
# ====== Computation of Mean  ===== #
# ==== Given Intermediate Fits ==== #
# ================================= #
Rcpp::cppFunction('
  Rcpp::NumericVector compute_mu_vec(const Rcpp::NumericVector& y,
                                     const Rcpp::NumericVector& y_tilde,
                                     const Rcpp::NumericVector& phi,
                                     const Rcpp::NumericVector& theta) {

    const int   n = y.size();
    const int   p = phi.size();
    const int   q = theta.size();
    Rcpp::NumericVector mu(n), tau(n);

    for (int t = 0; t < n; ++t) {

      // ------- AR part: 
      double ar = 0.0;
      for (int j = 0; j < p; ++j) {
        int idx = t - j - 1;             
        if (idx >= 0)
          ar += phi[j] * (y[idx] - y_tilde[idx]);
      }

      // ------- MA part: 
      double ma = 0.0;
      for (int j = 0; j < q; ++j) {
        int idx = t - j - 1;
        if (idx >= 0)
          ma += theta[j] * (y[idx] - mu[idx])/(mu[idx]);
      }

      tau[t] = ar + ma;
      mu[t]  = y_tilde[t] + tau[t];
    }
    return mu;
  }')


# ================================= #
# ====== Estimation of Shape  ===== #
# ==== Given Response and Fits ==== #
# ================================= #
compute_shape <- function(yt, mu_t, params, warnings = FALSE){
  # get 'effective n' and p
  n <- sum(!is.na(yt) & !is.na(mu_t))
  p <- length(params) + 1
  # squared Pearson residuals and MoM shape estimate
  alpha <- (n - p) / sum( ((yt - mu_t) / (mu_t))^2, na.rm = TRUE )
  return(alpha)
}
# alternatively, fixed alpha
mle_alpha <- function(y) {
  fit <- fitdistr(y, "gamma")
  unname(fit$estimate["shape"])
}
# =========================== #
# ====== Log-Likelihood ===== #
# =========================== #
# objective function for optim 
error <- function(xt, yt, params, configs) {
  # predict response using object parameters
  preds <- predict_response_GKR(xt, yt, params, configs)
  # check for infinite or NA fits
  if (any(!is.finite(preds$mu_hat)) | any(is.na(preds$mu_hat))) {
    return(-1e10)
  } else {
    # check stationarity
    if (configs$p > 0){
      # get AR params
      phi <- params[(3*configs$k + 1):(3*configs$k + configs$p)] 
      # if any roots are within unit circle in C 
      if (!all(Mod( polyroot(c(1, -phi)) ) > 1)){
        return(-1e10)
      }
    }
    # check invertibility
    if (configs$q > 0){
      # get MA params
      theta <- params[(length(params) - configs$q):(length(params) - 1)] 
      # if any roots are within unit circle in C 
      if (!all(Mod( polyroot(c(1, theta)) ) > 1)){
        return(-1e10)
      }
    }
    alpha <- compute_shape(yt = yt, mu_t = preds$mu_hat, params = params)
    # otherwise return the log-likelihood using our density function
    sum(dmeangamma(
      yt = yt,
      shape = alpha, # alpha is always the last entry of Xi
      mu = preds$mu_hat,
      log = TRUE
    ))
  }
}

# ======================== #
# ====== Column-Wise ===== #
# = Convolution via FFT. = #
# ======================== #
convolve_columns_fft <- function(X, kernels) {
  X <- as.matrix(X)
  n <- nrow(X)
  m <- ncol(X)
  if (is.matrix(kernels)) {
    if (ncol(kernels) != m)
      stop("ncol(kernels) must equal ncol(X)")
    kernels <- lapply(seq_len(m), function(j)
      kernels[, j])
  }
  if (length(kernels) != m)
    stop("length(kernels) must equal ncol(X)")
  # basic checks
  bad <- vapply(kernels, function(k)
    any(is.na(k)), logical(1))
  if (any(bad))
    stop(sprintf("Kernel(s) contain NA: %s", paste(which(bad), collapse = ",")))
  k_len <- pmin(vapply(kernels, length, integer(1)), n)
  if (any(vapply(kernels, length, integer(1)) > n))
    warning("Some kernels exceed series length; trimmed to n for efficiency")
  
  kmax <- max(k_len)
  pad_len <- 2^ceiling(log2(n + kmax - 1L))
  
  # zero-pad X to pad_len rows.
  X_pad <- rbind(X, matrix(0, nrow = pad_len - n, ncol = m))
  
  # build a pad_len m kernel matrix with per-column kernels, zero-padded.
  K_mat <- matrix(0, nrow = pad_len, ncol = m)
  for (j in seq_len(m)) {
    kj <- kernels[[j]][seq_len(k_len[j])]
    K_mat[seq_along(kj), j] <- kj
  }
  
  # FFT per column, multiply in frequency domain, inverse FFT; take first n rows.
  FX <- mvfft(X_pad)
  FK <- mvfft(K_mat)
  Y  <- mvfft(FX * FK, inverse = TRUE) / pad_len
  S  <- Re(Y[seq_len(n), , drop = FALSE])
  dimnames(S) <- dimnames(X)
  S
}
# X <- matrix(c(2, 3, 1, 1, -1, 5), nrow = 3, ncol = 2)
# K <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, ncol = 2)
# S <- convolve_columns_fft(X, K)
# S

build_kernel <- function(logdelta, logsigma, type = "gamma") { 
  # transform back the kernel parameters to the original scale
  delta <- exp(logdelta)
  sigma <- exp(logsigma)
  if (type == "gamma") {
    # transform into shape and scale
    alpha <- delta # delta^2 / sigma^2
    theta <- sigma # sigma^2 / delta
    # approximate limits using and quantiles roughly matching +/- 3 sigma
    min_lag <- floor(qgamma(0.0015, shape = alpha, scale = theta))
    max_lag <- ceiling(qgamma(0.9985, shape = alpha, scale = theta))
    # define the integer sequence
    lags    <- seq.int(min_lag, max_lag)
    # compute discrete kernel
    kernel        <- diff(pgamma(lags, shape = alpha, scale = theta))
    kernel[is.infinite(kernel)] <- 1
  }else if (type == "triangle") {
    # Compute the lower and upper limit of non-zero lags
    min_lag <- max(floor(delta - 3 * sigma + 1), 0)
    max_lag <- ceiling(delta + 3 * sigma - 1)
    lags    <- seq.int(min_lag, max_lag)
    
    # Compute the triangle weights values for the kernel (not normalized)
    kernel  <- 3 * sigma - abs(lags - delta)
  } else if (type == "gaussian") { # TODO: have two sets of gaussian kernels
    #      change to be equivalent across log-scales (adjust delta range)
    # define kernel radius by r = 3*sigma 
    kernelRad     <- max(3*sigma, 1) 
    min_lag <- max(floor(delta - 3 * sigma), 0)
    max_lag <- ceiling(delta + 3 * sigma)
    # define minimum / maximum index
    iKernel       <- c(min(floor(delta - kernelRad), 0), 
                       ceiling(delta + kernelRad)) 
    bins          <- seq(iKernel[1] - 0.5, iKernel[2] + 0.5)  # define bins
    # compute kernel values for window by integrating over Gaussian distribution in each bin
    kernel        <- diff(pnorm(bins, mean = delta, sd = sigma))
  }

  # # Pre-allocate and assign kernel values
  # kernel_full <- numeric(max_lag + 1)
  # kernel_full[(min_lag + 1):(max_lag + 1)] <- kernel

  # normalize the kernel to ensure it sums to 1
  kernel /  sum(kernel, na.rm = TRUE) 
}
# =========================== #
# ====== Hamons' Method ===== #
# = Of Daily PET from Temp. = #
# =========================== #
temp_to_PET <- function(temp, calendar_date, lat, long, startdate = NULL){
  # get sunlight hours from input dates
  sun_times <- getSunlightTimes(date = calendar_date, lat = lat, lon = long)
  # from this, compute daylight hours (sunrise to sunset)
  D  <- as.numeric(difftime(sun_times$sunset, sun_times$sunrise, units = "hours"))
  D  <- D /12
  # assign temperature to a new object
  Td <- temp
  # compute vapour pressure of water at t (Tetens equation)
  es <- sapply(Td, function(t){
    if(t > 0){
      0.61078*exp(17.27*t / (t + 237.3))
    } else{
      0.61078*exp(21.875*t / (t + 265.5))
    }
  })
  # convert es to Pa (multiply by 1000)
  es <- es * 1000  # convert kPa to Pa for consistency with SI units
  R_v <- 461.5  # J/(kgK)
  # re-do some conversions and finalize Pt
  Pt <- es / (R_v * (Td + 273.15)) 
  # compute slope vapour pressure unit temp
  slope_es <- mapply(function(e_s_value, t) {
    # convert es back to kPa 
    e_s_kPa <- e_s_value / 1000  
    # piecewise by temperature
    if (t > 0) {
      slope <- e_s_kPa * (17.27 * 237.3) / ((t + 237.3)^2)
    } else {
      slope <- e_s_kPa * (21.875 * 265.5) / ((t + 265.5)^2)
    }
    return(slope)
  }, es, Td)
  # constant
  r <- 0.0105
  # switch Pt to the correct units
  Pt_gm_m3 <- Pt * 1000 
  # convert everything to imperial
  Td_F <- Td * (9 / 5) + 32
  d_inHg <- slope_es * 0.2953 / 1.8
  # compute PET by Equation 3
  Ep_imp <- (3/4) * ((0.01 * r * D * Pt_gm_m3 + ((d_inHg * 5 * (D * Td_F - 27)) / 1500)) / (d_inHg + r))
  # compute PET by Equation 4
  C <- 0.0055
  Ep_imp4 <- C*(D^2)*Pt_gm_m3 
  # back to mm
  Ep_3 <- Ep_imp * 25.4
  Ep_4 <- Ep_imp4 * 25.4
  return(list(eq3 = Ep_3, eq4 = Ep_4))
}
