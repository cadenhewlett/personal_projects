// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=1> n;          
  vector[n] x;             
  vector[n] y;             
}

parameters { 
  // hyper-prior
  real<lower=0> s_0;
  real<lower=0> s_1;
  real<lower=0> s_2;
  real<lower=0> s_3;
  
  // coefficients for mean
  real b_0;
  real b_1;

  // coefficients for log-variance
  real g_0;
  real g_1;
}


model {
  // --------------- //
  // - Hyperpriors - //
  // --------------- //
  s_0 ~ exponential(1);
  s_1 ~ exponential(1);
  s_2 ~ exponential(1);
  s_3 ~ exponential(1);
  // ------------ //
  // --- Mean --- //
  // ------------ //
  b_0 ~ normal(0, s_0);
  b_1 ~ normal(0, s_1);
  // ---------------- //
  // - Log-Variance - //
  // ---------------- //
  g_0 ~ normal(0, s_2);
  g_1 ~ normal(0, s_3);
  // -------------- //
  // - Likelihood - //
  // -------------- //
  for (i in 1:n) {
    // mean
    real mu_i = b_0 + b_1 * x[i];
    // sigma = exp( (g_0 + g_1*x[i]) / 2 )
    real sigma_i = exp(0.5 * (g_0 + g_1 * x[i]));
    // normal likelihood
    y[i] ~ normal(mu_i, sigma_i);
  }
}


