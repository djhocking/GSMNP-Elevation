// Binomial mixture model with covariates
data {
  int<lower=0> R;       // Number of transects
  int<lower=0> nsites;           // Number of sites
  int<lower=0> T;       // Number of temporal replications
  int<lower=0> y[R, T]; // Counts
  vector[R] X;          // Covariate
  matrix[R, T] RH;      // Covariate
  int<lower=0> K;       // Upper bound of population size
  int<lower=0> sites[R];   // Sites 
  
}

transformed data {
  int<lower=0> max_y[R];

  for (i in 1:R)
    max_y[i] = max(y[i]);
}

parameters {
  real alpha0;
  real alpha1;
  real beta0;
  real beta1;
  vector[nsites] eps;            // Random site effects
  real<lower=0,upper=1> sd_eps;
  matrix[nsites, T] delta;            // Random site effects
  real<lower=0,upper=1> sd_p;
}

transformed parameters {
  vector[R] log_lambda; // Log population size
  matrix[R, T] logit_p; // Logit detection probability

for (i in 1:R) {
  log_lambda[i] = alpha0 + alpha1 * X[i] + eps[sites[i]];
  for (t in 1:T) {
  logit_p[i,t] = beta0 + beta1 * RH[i,t] + delta[sites[i], t];
  }
}
}

model {
  // Priors
  // Improper flat priors are implicitly used on
  // alpha0, alpha1, beta0 and beta1.

  eps ~ normal(0, sd_eps);
  sd_eps ~ uniform(0, 3);
  //  sd_eps ~ uniform(0, 1);   // Implicitly defined
  
  eps ~ normal(0, sd_p);
  sd_p ~ uniform(0, 3);
  
  // Likelihood
  for (i in 1:R) {
    vector[K - max_y[i] + 1] lp;

    for (j in 1:(K - max_y[i] + 1))
      lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
             + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]);
    target += log_sum_exp(lp);
  }
}

generated quantities {
  int N[R];
  int totalN;

  for (i in 1:R)
    N[i] = poisson_log_rng(log_lambda[i]);
  totalN = sum(N);
}