// Binomial mixture model with covariates
data {
  int<lower=0> R;       // Number of transects
  int<lower=0> T;       // Number of temporal replications
  int<lower=0> y[R, T]; // Counts
  vector[R] elev;          // Covariate
  vector[R] elev2;          // Covariate
  vector[R] litter;          // Covariate
  matrix[R, T] RH;      // Covariate
  matrix[R, T] temp;      // Covariate
  matrix[R, T] temp2;      // Covariate
  vector[R] gcover;      // Covariate
  vector[R] gcover2;      // Covariate
  int<lower=0> K;       // Upper bound of population size
}

transformed data {
  int<lower=0> max_y[R];
  int<lower=0> N_ll;
  int foo[R];

  for (i in 1:R)
    max_y[i] = max(y[i]);
    
  for (i in 1:R) {
    foo[i] = K - max_y[i] + 1;
}
    N_ll = sum(foo);
}
    

parameters {
  real alpha0;
  real alpha1;
  real alpha2;
  real alpha3;
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
}

transformed parameters {
  vector[R] log_lambda; // Log population size
  matrix[R, T] logit_p; // Logit detection probability

  log_lambda = alpha0 + alpha1 * elev + alpha2 * elev2 + alpha3 * litter;
  logit_p = beta0 + beta1 * RH + beta2 * temp + beta3 * temp2 + rep_matrix(beta4 * gcover, T);
}

model {
  // Priors
  // Improper flat priors are implicitly used on
  // alpha0, alpha1, beta0 and beta1.

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
  int<lower=0> N[R];                // Abundance
  int totalN;
  vector[R] log_lik;
  real mean_abundance;
  real fit = 0;
  real fit_new = 0;
  matrix[R, T] p; 
    matrix[R, T] eval;         // Expected values
    int y_new[R, T];
    matrix[R, T] E;
    matrix[R, T] E_new;
    vector[K + 1] lp;

  
    for (i in 1:R) {
    vector[K - max_y[i] + 1] ll;

    for (j in 1:(K - max_y[i] + 1)) {
      ll[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
             + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]);
    }
    log_lik[i] = log_sum_exp(ll);
  }
     
  for (i in 1:R) {
     N[i] = poisson_log_rng(log_lambda[i]);
     p[i, 1:T] = inv_logit(logit_p[i, 1:T]);
  }
  
  // Bayesian p-value fit

    // Initialize N, E and E_new
    // N = 0;
    for (i in 1:1) {
      for(j in 1:T) {
      E[i, j] = 0;
    E_new[i, j] = 0;
      }
    }
    for (i in 2:R) {
      E[i] = E[i - 1];
      E_new[i] = E_new[i - 1];
    }

    for (i in 1:R) {
      // for (j in 1:T)
        // p[i, j] = inv_logit(logit_p[i, j]);

        for (n in 0:(max_y[i] - 1))
          lp[n + 1] = negative_infinity();
        for (n in max_y[i]:K) {
          lp[n + 1] = poisson_log_lpmf(n | log_lambda[i])
            + binomial_lpmf(y[i, 1:T] | n, p[i, 1:T]);
        }
    
      for (j in 1:T) {
      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
      eval[i, j] = p[i, j] * N[i];
      E[i, j] = square(y[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
      // Generate replicate data and
      // Compute fit statistic E_new for replicate data
          y_new[i, j] = binomial_rng(N[i], p[i, j]);
          E_new[i, j] = square(y_new[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
        }
      }
    
   totalN = sum(N);  // Total pop. size across all sites
    for (i in 1:R) {
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
  mean_abundance = exp(alpha0);
}