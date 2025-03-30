library(tidyverse)
library(lmtest)
require(tidybayes)
require(rstan)
library(here)
# ----------- #
# -- SETUP -- #
# ----------- #
set.seed(1928)
# df <- read.csv("C:/Users/thetr/OneDrive/Desktop/Work/TA_STAT_200/pokedex.csv")
# generate heteroskedastic data
rstan_options(auto_write = TRUE)
x <- rnorm(100, mean = 0, sd = 1)
y <- 3*x - 4 + rnorm(100, mean = 0, sd = 6 + 1.5*x)
# format as data frame
df <- tibble(x = x, y = y)
n <- length(y)
# run HMC
here::here()
mod_b <- stan(
  "mod_b.stan",
  model_name = "mod_b",
  seed = 1928,
  data =   list(n = n, x = x, y = y), #compose_data(df),
  verbose = TRUE, # FALSE = not reporting,
  include = TRUE,
  chains = 4, 
  iter = 10000,   
  warmup = 3000
)
?sampling
saveRDS(mod_b, "mod_b.rds")
