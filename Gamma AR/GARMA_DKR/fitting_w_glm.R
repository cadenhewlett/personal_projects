library(DKR)
library(glarma)
source(here::here("Gamma AR", "GARMA_DKR", "garma_dkr_main.R"))
# experiment grounds
train <- koksilah[koksilah$hydr_year < 10, ]
set.seed(1928)
xt <- matrix(
  data = c(train$rain, train$rain, train$rain, train$rain*(train$rain > quantile(train$rain, 0.95))),
  ncol = 4,
  byrow = FALSE
)
yt <- train$gauge
# rough configs
configs <- list(k = 4, kernel_type = c("gamma", "gamma", "gamma", "gamma"))

# ------------------- #
# ---- Inner GLM ---- #
# ------------------- #
error <- function(xt, yt, kpars, configs, betas = NULL) {
  # then extract relevant hyper-parameters
  k <- configs$k
  kpars <- as.data.frame(matrix(kpars, ncol = 2))
  colnames(kpars) <- c("logdelta", "logsigma")
  betas <- NULL
  # build kernels for each set of parameters
  kernels <- sapply(1:k, function(i) {
    build_kernel(kpars$logdelta[i], kpars$logsigma[i], type = configs$kernel_type[i])
  })
  # convert kernel to list if there is only one kernel
  if (k == 1) {
    kernels <- list(kernels)
  }
  # convolve xt with kernels
  xt_conv <- convolve_columns_fft(xt, kernels = kernels)
  # prep convolved xts for IRLS
  subdf <- data.frame(y = yt, xt_conv)
  ols   <- lm(y ~ xt_conv, data = subdf)
  mu0   <- pmax(1e-6, as.numeric(predict(ols)))
  result <- tryCatch(
    suppressWarnings(glm(
      y ~ xt_conv,
      data = subdf,
      mustart = mu0,
      start = abs(coefficients(ols)),
      family = Gamma(link = "identity")
    )),
    error = function(e) NULL
  )
  # get log-likelihood
  if (is.null(result)) {
    return(list(loglik = -1e10))
  } else{
    betas <- coefficients(result)
    alpha <- 1 / summary(result)$dispersion
    ll <- sum(log(
      dmeangamma(
        yt = subdf$y,
        mu = result$fitted.values,
        shape = alpha
      )
    ))
    return(
      list(loglik = ll,
           kernel_pars = kpars,
           betas = betas, 
           alpha = alpha)
    )
  }
}

library(rgenoud)
ranges <- matrix(data = c(rep(-3, times = 2*configs$k), 
                          rep( 3, times = 2*configs$k)),
                 nrow = 2*configs$k,
                 byrow = FALSE)


set.seed(1928)
test <- genoud(
  fn = function(par){
    -1*error(xt = xt, yt = yt, kpars = par, configs = configs)$loglik
  },
  nvars           = 2*configs$k,
  starting.values = build_inits(n_inits = 100*configs$k, ranges = (ranges)),
  Domains         = (ranges),
  print.level     = 1,
  pop.size        = 100*configs$k, 
  boundary.enforcement = 2,
  max.generations = 50,
  BFGS = FALSE
)
build_kernel(3,3)


# ---- what does it look like? ---- # 
out <- error(xt = xt, yt = yt, kpars = test$par, configs = configs)
kpars <- out$kernel_pars
betas <- out$betas
alpha <- out$alpha

kernels <- sapply(1:k, function(i) {
  build_kernel(kpars$logdelta[i], kpars$logsigma[i], type = configs$kernel_type[i])
})

# convert kernel to list if there is only one kernel
if (k == 1) {
  kernels <- list(kernels)
}
# convolve xt with kernels
xt_conv <- convolve_columns_fft(xt, kernels = kernels)
# prep convolved xts for IRLS
subdf <- data.frame(y = yt, xt_conv)
ols   <- lm(y ~ xt_conv, data = subdf)
mu0   <- pmax(1e-6, as.numeric(predict(ols)))
result <- tryCatch(
  suppressWarnings(glm(
    y ~ xt_conv,
    data = subdf,
    mustart = mu0,
    start = abs(coefficients(ols)),
    family = Gamma(link = "identity")
  )),
  error = function(e) NULL
)
# get log-



mu_pred <- result$fitted.values
y_obs <- yt
unit_dev <- 2  * ( (y_obs - mu_pred)/mu_pred - log(y_obs/mu_pred) )
dev_resd <- sign(y_obs - mu_pred)*sqrt(unit_dev)

par(mfrow = c(1,2))
qqnorm(dev_resd[!is.na(dev_resd)], main = "Deviance Residual Q-Q Plot")
abline(0, 1, col = 'red')
plot(result$fitted.values[250:950], type = 'l', ylim = range(yt),
     main = "Actual (Red) vs. Predicted (Black)", ylab  = "Streamflow")
lines(yt[250:950], col = 'red')

