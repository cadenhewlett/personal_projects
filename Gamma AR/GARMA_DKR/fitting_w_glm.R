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
xt <- matrix(
  data = c(rep(1, times = nrow(train)), 
           train$rain, 
           train$rain,
           train$temp_kelvin
           ),
  ncol = 4,
  byrow = FALSE
)
yt <- train$gauge

cap <- quantile(train$gauge, 0.75) + 1.5*IQR(train$gauge)
w <- 1 + train$gauge*(train$gauge >= cap)
# rough configs
configs <- list(k = ncol(xt),
                kernel_type = c("fixed", "gamma",  "gamma", "gamma"),
                thresholded = NULL,
                c = 0,
                q = 1,
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

predict_response_GKR <- function(xt, yt, pars, configs, betas){
  # then extract relevant hyper-parameters
  k <- configs$k 
  # p <- configs$p
  # q <- configs$q
  # extract number of fixed kernels
  n_fixed <- sum(configs$kernel_type == "fixed")
  # extract kernel params
  kpars <- rbind(matrix(
    rep(0, times = 2 * n_fixed), ncol = 2), 
    as.data.frame(matrix(pars[1:(2*(k-n_fixed))], ncol = 2))
  )
  colnames(kpars) <- c("logdelta", "logsigma")
  betas <- NULL
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
  # apply thresholding
  xt_conv[, configs$thresholded] <- sapply(1:configs$c, function(i) {
    idx <- configs$thresholded[i]
    xt_conv[, idx] * (xt_conv[, idx] >  pars[2 * (k - n_fixed) + i])
  })
  if(is.null(betas)){
    betas <- irls_gamma_identity_qr(X = xt_conv, y = yt, w = configs$w)
  }
  mu_fits <- drop(xt_conv %*% betas)
  # compute alpha
  alpha <- compute_shape(yt = yt, mu_t = mu_fits, params = c(betas, pars))
  # get loglik
  ll <- sum(log(
    dmeangamma(
      yt = yt,
      mu = mu_fits,
      shape = alpha
    )
  ))
  return(
    list(loglik = ll,
         kernel_pars = kpars,
         thresholds =  pars[2 * (k - n_fixed) + 1 : 2 * (k - n_fixed) + configs$c],
         betas = betas, 
         alpha = alpha,
         mu_fits = mu_fits)
  )
}

error <- function(xt, yt, pars, configs, betas = NULL) {
  # then extract relevant hyper-parameters
  k <- configs$k 
  p <- configs$p
  q <- configs$q
  c <- configs$c
  n_fixed <- sum(configs$kernel_type == "fixed")
  # extract kernel params
  kpars <- rbind(matrix(
    rep(0, times = 2 * n_fixed), ncol = 2), 
    as.data.frame(matrix(pars[1:(2*(k-n_fixed))], ncol = 2))
  )
  colnames(kpars) <- c("logdelta", "logsigma")
  betas <- NULL
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
  # apply thresholding
  if(configs$c > 0){
    xt_conv[, configs$thresholded] <- sapply(1:configs$c, function(i) {
      idx <- configs$thresholded[i]
      xt_conv[, idx] * (xt_conv[, idx] >  pars[2 * (k - n_fixed) + i])
    })
  }
  # use RCpp QR procedure to estimate betas
  betas <- irls_gamma_identity_qr(X = xt_conv, y = yt, w = configs$w)
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
  ll <- sum(log(
      dmeangamma(
        yt = yt,
        mu = mu_t,
        shape = alpha
      )
    ))
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


library(rgenoud)

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
  rep(c(-0.5, 0.5), times = configs$q),
  # alpha
  c(0, 20)
)

set.seed(1928)
fit <- genoud(
  fn = function(par){
    -1*error(xt = xt, yt = yt, pars = par, configs = configs)$loglik
  },
  nvars           = nrow(ranges),
  starting.values = build_inits(n_inits = 100*nrow(ranges), ranges = ranges),
  Domains         = ranges,
  print.level     = 1,
  pop.size        = 100*nrow(ranges), 
  boundary.enforcement = 2,
  max.generations = 20,
  BFGS = FALSE
)
fit$par
test_set <- data[data$hydr_year %in% c(11, 12), ]
set.seed(1928)
out <- error(xt, yt, fit$par, configs = configs)
xt <- matrix(
  data = c(rep(1, times = nrow(test_set)),
           test_set$rain,
           test_set$rain,
           test_set$temp_kelvin),
  ncol = 4,
  byrow = FALSE
)
yt <- test_set$gauge
# test_set_configs <- configs
# test_set_configs$w <- 1 + (test_set$gauge >= cap)

kpars <- out$kernel_pars
betas <- out$betas
# then extract relevant hyper-parameters
k <- configs$k 
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
n_fixed <- sum(configs$kernel_type == "fixed")
# apply thresholding
# xt_conv[, configs$thresholded] <- sapply(1:configs$c, function(i) {
#   idx <- configs$thresholded[i]
#   xt_conv[, idx] * (xt_conv[, idx] >  fit$par[2 * (k - n_fixed) + i])
# })
# prep convolved xts for IRLS
mu_comp <- drop( xt_conv %*% out$betas )
mu_pred <- compute_mu_vec(y = yt, y_tilde = mu_comp, phi = fit$par[7], theta = fit$par[8])
# get log-
plot(mu_pred, type = 'l')
lines(yt, col = 'red')
y_obs <- yt
unit_dev <- 2  * ( (y_obs - mu_comp)/mu_comp - log(y_obs/mu_comp) )
dev_resd <- sign(y_obs - mu_comp)*sqrt(unit_dev)

pearson <- ((y_obs - mu_comp)/sqrt(mu_comp^2/fit$par[length(fit$par)])) 
qqnorm(pearson)
abline(0, 1, col = 'red')
par(mfrow = c(1,2))
qqnorm(dev_resd[!is.na(dev_resd)], main = "Deviance Residual Q-Q Plot")

abline(0, 1, col = 'red')
plot(mu_pred[150:750], type = 'l', ylim = range(yt),
     main = "Actual (Red) vs. Predicted (Black)", ylab  = "Streamflow")

lines(yt[150:750], col = 'red')
# r2
r2 <- 1-sum( (yt - mu_pred)^2)/sum ( (yt - mean(yt))^2)
kge <- hydroGOF::KGE(sim = mu_pred, obs = yt)
nse <- hydroGOF::NSE(sim = mu_pred, obs = yt)
round(c(r2, nse, kge), digits = 2)

