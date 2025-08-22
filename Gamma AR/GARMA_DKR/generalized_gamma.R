# GENERALIZED GAMMA VERSION
predict_response_GKR <- function(xt, yt, params, configs) {
  # default config is no AR, no MA
  configs_default <- list(
    k = max(1, length(params)/3 - 1), 
    p = 0,  q = 0, g = I, c = 0,
    kernel_type = "gamma",
    thresholded = NULL,
    use_log = FALSE 
  )
  # modify by user configuration
  configs <- modifyList(configs_default, configs)
  # then extract relevant hyper-parameters
  k <- configs$k
  p <- configs$p
  q <- configs$q
  # ---------- #
  # --- KR --- #
  # ---------- #
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
  theta <- if (q > 0) {
    params[(length(params) - 2 - q + 1):(length(params) - 2)]
  } else numeric(0)
  # compute GARMA adjusted mean
  mu_t <- compute_mu_vec(y = yt, y_tilde = yt_tilde, phi, theta)
  if (any(!is.finite(mu_t)) || any(mu_t <= 0)) {
    return(list(loglik = -1e10, mu_hat = mu_t))
  }
  #  family parameters
  kappa <- exp(params[length(params) - 1])
  tau   <- exp(params[length(params)])  
  # compute loglik
  ll <- sum(dggamma_mean(yt, mu = mu_t, kappa = kappa, tau = tau, log = TRUE))
  # return results
  list(
    loglik         = ll,
    mu_hat         = mu_t,
    yt_tilde       = yt_tilde,
    tau_t          = mu_t - yt_tilde,
    kappa          = kappa,
    tau            = tau,
    kpars          = pars[, c("logdelta", "logsigma")],
    beta           = pars$beta
  )
}

# generalized gamma density with mean parameterization
dggamma_mean <- function(y, mu, kappa, tau, log = FALSE) {
  # enforce domain; return -Inf if invalid to help optimizers
  if (any(mu <= 0) || kappa <= 0 || tau <= 0) {
    out <- rep(-1e10, length(y)); return(if (log) out else exp(out))
  }
  # scale implied by mean
  mratio <- gamma((kappa + 1)/tau) / gamma(kappa/tau)   # m(kappa,tau)
  beta   <- mu / mratio
  # log pdf (stable)
  logf <- log(tau) - lgamma(kappa/tau) - kappa*log(beta) +
    (kappa - 1)*log(y) - (y/beta)^tau
  logf[!(y > 0)] <- -Inf
  if (log) logf else exp(logf)
}

# ------------------- #
# -- Configuration -- #
# ------------------- #
# info on koksilah river from google maps
lat <- 48.75603116961046
long <-  -123.6820397210659
data <- DKR::koksilah
# configure dates
start_date <- as.Date("01-10-1988", format = "%d-%m-%Y")
days_offset <- seq_along(data$hydr_year) - 1
calendar_date <- start_date + days_offset
# compute transform and append
data$PET <- temp_to_PET(data$temp, calendar_date, lat, long)$eq3
data$temp_kelvin <- data$temp + 273.15
# experiment grounds
train <- data[data$hydr_year < 30, ]
set.seed(1928)
covariates <- c(rep(1, times = nrow(train)), 
                train$rain, 
                train$rain,
                train$rain,
                train$rain*train$temp_kelvin
)
xt <- matrix(
  data = covariates,
  ncol = length(covariates) / nrow(train),
  byrow = FALSE
)
yt <- train$gauge
configs$k <- ncol(xt)
configs$p <- 0
configs$q <- 0
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
    -3, 3,
    -3, 3
  ),
  ncol = 2,
  byrow = TRUE)
)
ranges <- R

fit <- genoud(
  # fn = function(par){
  #   -1*predict_response_GKR(xt = xt, yt = yt, params = par, configs = configs)$loglik
  # },
  fn = function(par){
    lambda <- 0.25
    fit <- predict_response_GKR(xt, yt, par, configs)
    if (!is.finite(fit$loglik)) return(1e12)
    # mixed objective: maximize likelihood AND minimize SSE
    sse <- sum((yt - fit$mu_hat)^2)
    -( fit$loglik - lambda * sse )         # choose lambda by CV
  },
  nvars           = nrow(ranges),
  starting.values = build_inits(n_inits = 100*nrow(ranges), ranges = ranges),
  Domains         = ranges,
  print.level     = 1,
  pop.size        = 100*nrow(ranges), 
  boundary.enforcement = 2,
  max.generations = 50,
  BFGS = FALSE
)

mu_pred <- predict_response_GKR(xt, yt, params = fit$par, configs)$mu_hat
y_obs <- yt

hydroGOF::KGE(sim = mu_pred, obs = y_obs)

fit$par

#  family parameters
kappa <- exp(fit$par[length(fit$par) - 1])
tau   <- exp(fit$par[length(fit$par)])  
# CDF for GG with mean parameterization
pggamma_mean <- function(y, mu, kappa, tau){
  m1 <- gamma((kappa+1)/tau) / gamma(kappa/tau)
  beta <- mu / m1
  pgamma( (y/beta)^tau, shape = kappa/tau, scale = 1 )
}



# test set 
test <- data[data$hydr_year >= 30, ]
covariates <- c(rep(1, times = nrow(test)), 
                test$rain, 
                test$rain,
                test$rain,
                test$rain*test$temp_kelvin
)
xt_test <- matrix(
  data = covariates,
  ncol = length(covariates) / nrow(test),
  byrow = FALSE
)
yt_test <- test$gauge
out_test <- predict_response_GKR(xt = xt_test, yt = yt_test, params = fit$par, configs = configs )

plot(out_test$mu_hat, type = 'l'); lines(yt_test, col = 'red')

hydroGOF::NSE(sim = out_test$mu_hat, obs = yt_test)

u  <- pggamma_mean(yt_test, out_test$mu_hat, kappa, tau)
rQ <- qnorm(pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps))
par(mfrow=c(2,2))
hist(rQ, breaks=40, main="Dunn–Smyth residuals", xlab="")
qqnorm(rQ); abline(0,1, col = 'red')
plot(y_obs[1:550], type = 'l', main = "Actual (Black) vs. Fitted Mean (Red)"); lines(mu_pred[1:550], col = 'red')
plot(rQ ~ out_test$mu_hat, cex=.6, main="rQ vs fitted mean");abline(h = 0, col = 'red')

# slope & intercept of y on mu (ideal: a=0, b=1)
# cal <- lm(yt ~ 0 + mu_pred)           # or lm(yt ~ mu_pred)
# summary(cal)$coef
# brks <- quantile(mu_pred, probs = seq(0,1,0.1), na.rm=TRUE)
# bin  <- cut(mu_pred, brks, include.lowest=TRUE)
# calib <- aggregate(cbind(y=yt, mu=mu_pred) ~ bin, FUN=mean)
# plot(calib$mu, calib$y, pch=16); abline(0,1,col=2,lwd=2)
# # check ccf
# cc <- ccf(yt, mu_pred, lag.max = 20, plot=T)
# best_lag <- cc$lag[ which.max(abs(cc$acf)) ]
# best_lag 
# 
# # binwise calibration curve
# brks <- quantile(mu_pred, probs = seq(0,1,0.1), na.rm=TRUE)
# bin  <- cut(mu_pred, brks, include.lowest=TRUE)
# calib <- aggregate(cbind(y=yt, mu=mu_pred) ~ bin, FUN=mean)
# plot(calib$mu, calib$y, pch=16); abline(0,1,col=2,lwd=2)
