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
train <- data[data$hydr_year < 30, ]
# data matrix
xt <- matrix(
  data = c(
    train$rain , train$rain, train$rain*train$PET # train$rain*train$PET,
  ), ncol = 3, byrow = FALSE
)
yt <- train$gauge
# configurations
configs <- list(
  k = 3, 
  p = 1,  q = 0, g = 0,
  kernel_type = rep("gamma", times = 3),
  use_log = FALSE 
)
# fit a model
set.seed(1928)
fit <- fit_GKR(
  xt = xt, yt = yt, configs = configs, report = TRUE
)
# ------------------------ #
# --- Test Performance --- #
# ------------------------ #
test <- data[data$hydr_year >= 30, ]
# data matrix
test_xt <- matrix(
  data = c(
    train$rain , train$rain, train$rain*train$PET 
  ), ncol = 3, byrow = FALSE
)
test_yt <- test$gauge[1:which(is.na(test$gauge))[1]]
oos <- predict_response_GKR(xt = test_xt, yt = test_yt, 
                            params = fit$par, configs = configs)
# how does it perform
fitted <- oos$mu_hat[max(configs$p, configs$q):which(is.na(test_yt))[1]]
actual <-  test_yt[1:(which(is.na(test_yt))[1] - max(configs$p, configs$q) + 1)]
plot(oos$mu_hat, type = 'l')# look at residuals
par(mfrow = c(1, 2))
y <- actual[-1]
mu <- fitted[-1]
shape <- fit$par[length(fit$par)]
dev <- sign(y-mu) * sqrt(2*shape*(y/mu - log(y/mu) - 1))
pearson <- (y - mu) / sqrt(mu^2 / shape)

# 
qres <- qnorm(pgamma(y, shape = shape, scale = mu/shape))
qres <- qres[-which(is.na(qres))]
qres <- qres[-which(!is.finite(qres))]
qqnorm(qres); qqline(qres, col = 2, lwd = 2)

plot(mu, pearson,
     xlab = "Fitted mean", ylab = "Pearson residual",
     main = "Pearson Residuals vs. Fitted Mean");  abline(h = 0, col = 2)

sum( pearson >= 2, na.rm = T ) / length(pearson)
sum( pearson <= -2, na.rm = T ) / length(pearson)

residuals <- (y - mu)
Box.test(qres,  lag = 20, type = "Ljung")   # remaining autocorrelation
acf(qres)

plot(x = fitted, y =  actual - fitted)
abline(h = 0, col = 'red')
hydroGOF::KGE(obs = actual, sim = fitted)
hydroGOF::NSE(obs = actual, sim = fitted)

# notes
# ok... so the primary issue I have now is that the (original scale) metrics are laughably bad...
# 
# I am worried that there is no way at all to marry the idea of theoretical well-suited 
# residuals and practically useful fitted values and coefs for this particular data set. 