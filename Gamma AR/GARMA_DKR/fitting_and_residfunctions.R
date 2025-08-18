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
# -------------------- #
# ----- Training ----- #
# -------------------- #
# configure wrt to our example
train <- data[data$hydr_year < 15, ]
# data matrix
xt <- matrix(
  data = c(
    rep(1, times = nrow(train)), train$rain, train$rain, train$PET
  ), nrow = length(train$rain), byrow = FALSE
)
yt <- train$gauge
# test set
test <- data[data$hydr_year %in% 15:20, ]
test_xt <- matrix(
  data = c(
    rep(1, times = nrow(test)), test$rain, test$rain, test$PET
  ), nrow = length(test$rain), byrow = FALSE
)
test_yt <- test$gauge
# configurations
configs <- list(
  k = ncol(xt), 
  p = 0,  q = 0, g = 0,
  kernel_type = rep("gamma", times = ncol(xt)),
  use_log = FALSE
)
# fit a model
set.seed(1928)
fit <- fit_GKR(
  xt = xt, yt = yt, configs = configs, report = TRUE
)

# --------------------- #
# -- Compute Kernels -- #
# --------------------- #
kernels <- matrix(fit$par[1:(3 * fit$configs$k)], nrow = fit$configs$k, byrow = TRUE)
colnames(kernels) <- c("beta", "logdelta", "logsigma")

# -------------------------- #
# -- Compute Fitted Vals -- #
# -------------------------- #


fits <- predict_response_GKR(
  xt = test_xt,
  yt = test_yt,
  params = fit$par,
  configs = fit$configs
)
dmeangamma(yt = test_yt, shape = fit$alpha_hat, mu = fits$mu_hat, log = TRUE)
ell    <- max(configs$p, configs$q)
n <- length(test_yt)
mu_pred   <- fits$mu_hat[(ell+1):n]
y_obs <- test_yt[(ell+1):n]

plot(mu_pred, type = 'l', col = 'red', ylim = range(yt))
lines(y_obs)
# --------------------------------------- #
# -- Residual Diagnostics: Computation -- #
# --------------------------------------- #
n        <- length(y_obs)
p        <- length(fit$par) + 1
alpha    <- fit$alpha_hat

alpha <- compute_shape(y_obs, mu_pred, c(fit$par, 1))
pearson  <- sqrt(alpha) * (y_obs - mu_pred) / mu_pred
# tinkering around
qqnorm(1.*pearson)
abline(0, 1)
# GoF statistic
# deviance residuals
unit_dev <- 2  * ( (y_obs - mu_pred)/mu_pred - log(y_obs/mu_pred) )
dev_resd <- sign(y_obs - mu_pred)*sqrt(unit_dev)

par(mfrow = c(1,2))
qqnorm(dev_resd[!is.na(dev_resd)], main = "Deviance Residual Q-Q Plot")
abline(0, 1)


FD      <- (fit$alpha_hat / (n - p))*sum(unit_dev)
FP      <- fit$alpha_hat * sum(pearson^2) / (n - p)
c(fit$par[-(1:(3*configs$k))], fit$alpha_hat)

# F_{\,\text{D}} =
# \frac{\hat{\alpha}}{(n-p)}\sum_{i = 1}^n d(y_i, \hat{\mu}_i)  \sim F_{(n-p), \infty}

kernels
fit$par[-(1:(3*configs$k))]


hydroGOF::NSE(sim = mu_pred, obs = y_obs)
hydroGOF::KGE(sim = mu_pred, obs = y_obs)
