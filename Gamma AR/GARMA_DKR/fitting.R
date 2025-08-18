source(here::here("Gamma AR", "GARMA_DKR", "garma_dkr_main.R"))
# get data
set.seed(1928)
library(DKR)
library(pso)
data <- DKR::koksilah_PET
train <- data[data$hydr_year < 30,]
configs <- list(
  k=2,p=2,q=2
)
xt <- train$rain
zt <- train$PET
yt <- train$gauge
garma_fit <- fit_DKR(
  xt = xt,
  zt = zt,
  yt = yt,
  configs = configs,
  report = TRUE
)


test <- data[data$hydr_year >= 30, ]

# -------------- #
# -- Evaluate -- #
# -------------- #
p <- configs$p
q <- configs$q
preds <- garma_fit$mu_t
plot(preds[(2):150], type = 'l')
lines(yt[1:(149)], col = 'red')
# compute resids
residuals <- yt[1:(length(yt)-1)]-preds[2:length(yt)]
DKR::eval_all(yt = yt[1:(length(yt)-1)], 
              yt_hat = preds[2:length(yt)])

plot( y = residuals[200:450], preds[200:450])

# ----------------- #
library(ggplot2)
library(dplyr)
library(tidyr)

# prepare data
df <- data.frame(
  time = 365:730,
  tau = garma_fit$tau_t[365:730],
  ytilde = garma_fit$yt_tilde[365:730]
) %>%
  mutate(mu = tau + ytilde)

# long format for tau and ytilde (for stacking)
df_long <- df %>%
  select(time, tau, ytilde) %>%
  pivot_longer(cols = c(tau, ytilde),
               names_to = "component",
               values_to = "value") %>%
  mutate(component = factor(component, levels = c("tau", "ytilde")))  
# plot
ggplot() +
  geom_area(
    data = df_long,
    aes(x = time, y = value, fill = component),
    position = "stack",
    alpha = 0.5
  ) +
  geom_line(
    data = df,
    aes(x = time, y = mu),
    color = "black",
    size = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    title = expression("Decomposition of " * mu[t] * " into " * tau[t] * " and " * tilde(y)[t]),
    x = "Time",
    y = "Value",
    fill = "Component"
  ) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c(expression(tau[t]), expression(tilde(y)[t]))) +
  theme_bw() + theme(legend.position = "bottom")


# ------------------- #
# -- Out of Sample -- #
# ------------------- #,
par(mfrow = c(1,1))
test <- data[data$hydr_year >= 30,]
test_xt <- test$rain
test_yt <- test$gauge
test_zt <- test$PET

# fit oos
test_fit <- predict_response_garma(test_xt, test_zt, test_yt, garma_fit$par, configs)

# a bit *too* good..?
# predicted vs. fitted
plot(test_fit$mu_hat[max(p,q):1000], type = 'l',
     main = "Predicted vs. Actual")
lines(test_yt[1:1000], col = 'red')

plot(y=test_yt[1:1000]-test_fit$mu_hat[max(p,q):1001], 
     x=test_fit$mu_hat[max(p,q):1001], type = 'p',
     main = "Residual vs. Fitted")
round(DKR::eval_all(yt = test_yt[1:1000], 
              yt_hat = test_fit$mu_hat[max(p,q):(1000+1)]),2)
  
