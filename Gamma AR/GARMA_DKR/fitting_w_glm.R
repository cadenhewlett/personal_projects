# load the helper functions
source(here::here("Gamma AR", "GARMA_DKR", "garma_dkr_main.R"))
# ------------------- #
# -- PET Transform -- #
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
                train$temp_kelvin
)
xt <- matrix(
  data = covariates,
  ncol = length(covariates) / nrow(train),
  byrow = FALSE
)

yt <- train$gauge
cap <- quantile(train$gauge, 0.75) + 1.5*IQR(train$gauge) 
cap <- quantile(train$gauge, 0.95)

# some weights for high yt 
w <- 1 + train$gauge*(train$gauge >= cap)
w <- w # / mean(w)
# rough configs
configs <- list(k = ncol(xt),
                kernel_type = c("fixed", "gamma",  "gamma", "gamma", "gamma"),
                thresholded = c(2),
                c = 1,
                q = 1  ,
                p = 1,
                w = w
                )


# ------------------- #
# ---- Inner GLM ---- #
# ------------------- #
# pars is of the form...
# [ kernel pars, threshold pars, ar pars, ma pars, shape ]
# let k = (number of paramaterized kernels) 
# pars[1:(2*k)] -> kernel params
# let c = (number of threshold variables)
# pars[(2*k)+1:(2*k)+c] -> threshold pars
# then pars[(2*k)+c+1 : (2*k)+c+p ] -> ar pars
# and (length(params) - q):(length(params) - 1) -> ma pars
# and pars[length(pars)] -> alpha

predict_response_GKR <- function(xt, yt, pars, configs, betas = NULL){
  # then extract relevant hyper-parameters
  k <- configs$k 
  p <- configs$p
  q <- configs$q
  c <- configs$c
  # compute max lag index
  ell <- max(p, q)
  # extract number of fixed kernels
  n_fixed <- sum(configs$kernel_type == "fixed")
  # extract kernel params
  kpars <- rbind(matrix(
    rep(0, times = 2 * n_fixed), ncol = 2), 
    as.data.frame(matrix(pars[1:(2*(k-n_fixed))], ncol = 2))
  )
  colnames(kpars) <- c("logdelta", "logsigma")
  # build kernels for each set of parameters
  kernels <- sapply(1:k, function(i) {
    if(configs$kernel_type[i] == "fixed"){
      build_kernel(0, 0, type = "fixed")
    } else{
      build_kernel(kpars$logdelta[i], kpars$logsigma[i], type = configs$kernel_type[i])
    }
  })
  # convert kernel to list if there is only one kernel
  if (k == 1) {
    kernels <- list(kernels)
  }
  # convolve xt with kernels
  xt_conv <- convolve_columns_fft(xt, kernels = kernels)
  # if no betas are provided
  if(is.null(betas)){
    # use RCpp QR procedure to estimate betas
    betas <- irls_gamma_identity_qr(X = xt_conv, y = yt, w = configs$w)
  }
  # then take prod
  mu_fits <- drop(xt_conv %*% betas)
  # extract phi if present
  phi <- if (p > 0) {
    pars[(2 * (k - n_fixed) + c + 1):(2 * (k - n_fixed) + c + p)]
  } else{
    numeric(0)
  }
  # ibid for theta
  theta <- if (q > 0) {
    pars[(length(pars) - q):(length(pars) - 1)]
  } else{
    numeric(0)
  }
  # GARMA errors
  mu_t <- compute_mu_vec(y = yt, y_tilde = mu_fits, phi, theta)
  # compute alpha
  alpha <- pars[length(pars)] #compute_shape(yt = yt, mu_t = mu_fits, params = c(betas, pars))
  # get loglik
  ll <- sum((dmeangamma(
    yt = yt,
    mu = mu_t,
    shape = alpha,
    log = TRUE
  )))
  return(
    list(loglik = ll,
         kernel_pars = kpars,
         pars = pars,
         betas = betas, 
         alpha = alpha,
         mu_fits = mu_t,
         mu_raw  = mu_fits)
  )
}

error <- function(xt, yt, pars, configs, betas = NULL) {
  # use fitting function and get log-likelihood
  predict_response_GKR(xt = xt, yt = yt, pars = pars, 
                       configs = configs, betas = NULL)$loglik
  
}


library(rgenoud)

fit_GKR <- function(xt, yt, configs){
  
  configs$nkernels <- configs$k - sum(configs$kernel_type == "fixed")
  
  ranges <- rbind(
    # logsigma, logdelta ranges for free kernels
    do.call(rbind, lapply(1:configs$nkernels, function(i) {
      matrix(c(
        # lower
        -3, -3,
        # upper
        3, 3
      ), ncol = 2, byrow = FALSE)
    })),
    # thresholding in response (if present)
    (if (configs$c > 0) {
      matrix(c(
        sapply(configs$thresholded, function(c) {
          range(xt[, c])
        })[1, ],
        sapply(configs$thresholded, function(c) {
          range(xt[, c])
        })[2, ]
      ), ncol = 2, byrow = FALSE)
    } else{
      numeric(0)
    }),
    # ARMA
    rep(c(-0.5, 0.5), times = configs$p),
    rep(c(-0.75, 0.75), times = configs$q),
    # alpha
    c(0, 20)
  ) 
  
  # --------- fit 1: core model
  fit <- genoud(
    fn = function(par){
      -1*error(xt = xt, yt = yt, pars = par, configs = configs)
    },
    nvars           = nrow(ranges),
    starting.values = build_inits(n_inits = 100*nrow(ranges), ranges = ranges),
    Domains         = ranges,
    print.level     = 1,
    pop.size        = 100*nrow(ranges), 
    boundary.enforcement = 2,
    max.generations = 15,
    BFGS = FALSE
  )

  return(fit)
}


# --------- TESTING ZONE ------------ #
fits <- fit_GKR(xt, yt, configs = configs)
# fitted
fits$par
betas <- predict_response_GKR(xt = xt,
                            yt = yt,
                            pars = fits$par,
                            configs = configs, 
                            betas = NULL)$betas
# ---------- out of sample
test_set <- data[data$hydr_year >= 30, ]

xt_test <- matrix(
  data = c(rep(1, times = nrow(test_set)), 
           test_set$rain, 
           test_set$rain,
           test_set$rain,
           test_set$temp_kelvin),
  ncol = ncol(xt),
  byrow = FALSE
)
yt_test <- test_set$gauge[1:nrow(xt_test)]

idx <- min(which(is.na(yt_test)))-1
yt_test <- yt_test[1:idx]
xt_test <- xt_test[1:idx, ]

out <- predict_response_GKR(
  xt = xt_test,
  yt = yt_test,
  pars = fits$par,
  configs = configs,
  betas = betas
)



# diagnostics
y_obs <- yt_test
mu_pred <- pmax(0, out$mu_fits[1:idx])
plot(mu_pred, type = 'l')
lines(y_obs, col = 'red')

unit_dev <- 2  * ( (y_obs - mu_pred)/mu_pred - log(y_obs/mu_pred) )
dev_resd <- sign(y_obs - mu_pred)*sqrt(unit_dev)

par(mfrow = c(1,2))
qqnorm(dev_resd[!is.na(dev_resd)], main = "Deviance Residual Q-Q Plot")
abline(0, 1)

configs$nkernel <- configs$k - sum(configs$kernel_type == "fixed")

# 
acf(y_obs - mu_pred)
#
r2 <- 1-sum( (y_obs - mu_pred)^2 ) / sum( (y_obs - mean(mu_pred))^2 )
nse <- hydroGOF::NSE(sim = mu_pred, obs = y_obs)
kge <- hydroGOF::KGE(sim = mu_pred, obs = y_obs)
# 
print(round(c(r2, nse, kge), 3))


n <- length(y_obs)
alpha <- fit$par[length(fit$par)]

# deviance residuals


## ---------- (Randomized) Quantile / PIT residuals ----------

# For continuous Y, randomized and non-randomized are identical:
pit <- pgamma(y_obs, shape = alpha, scale = mu_pred/alpha)
eps <- .Machine$double.eps
pit <- pmin(pmax(pit, eps), 1 - eps)
# 
rqr <- qnorm(pit)  # quantile residuals
# Q–Q plot against N(0,1)
qqnorm(rqr, main = "Quantile residuals: Q–Q vs N(0,1)")
abline(0, 1, col = 'red')
# ACF of quantile residuals (should be white if model is calibrated)
acf(rqr, main = "ACF of quantile residuals")

# PIT histogram (should be ~flat if predictive distribution is calibrated)
hist(pit, breaks = "FD", freq = TRUE,
     main = "PIT histogram", xlab = "u = F(y_t | F_{t-1})")

## ---------- Standardized innovation residuals ----------
innov <- (y_obs - mu_pred) / sqrt(
  (mu_pred^2) / alpha)

# Quick moment checks
cat(sprintf("mean(innov)=%.3f, var(innov)=%.3f\n",
            mean(innov, na.rm=TRUE), var(innov, na.rm=TRUE)))

# ACF (raw) and ACF of squares (for leftover conditional heteroskedasticity)
acf(innov,     main = "ACF of standardized innovations")
acf(innov^2,   main = "ACF of squared standardized innovations")

# Ljung–Box on several lags (adjust 'lag' as appropriate to series length)
# H0: no autocorrelation up to 'lag'
for (L in c(10, 20, 30)) {
  lb <- Box.test(innov, lag = L, type = "Ljung-Box", fitdf = 0)
  cat(sprintf("Ljung–Box lag=%d: p=%.3g\n", L, lb$p.value))
}
