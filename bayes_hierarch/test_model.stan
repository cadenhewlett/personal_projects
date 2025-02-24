data {
  int<lower=0> N;          // Number of observations
  vector[N] y;             // Observed data
}

parameters {
  real mu;                 // Mean parameter
  real<lower=0> sigma;     // Standard deviation
}

model {
  mu ~ normal(0, 10);      // Prior for mean
  sigma ~ normal(0, 5);    // Prior for sigma
  y ~ normal(mu, sigma);   // Likelihood
}
