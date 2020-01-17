// Binomial mixture model with covariates
data {
  int<lower=0> R;       // Number of transects
  int<lower=0> T;       // Number of temporal replications
  int<lower=0> nsites;       // Number of sites
  int<lower=1> sites[R];       // vector of sites
  int<lower=0> y[R, T]; // Counts
  vector[R] elev;          // Covariate
  vector[R] elev2;          // Covariate
  vector[R] litter;          // Covariate
  vector[R] twi;          // Covariate
  vector[R] stream;          // Covariate
  matrix[R, T] RH;      // Covariate
  matrix[R, T] precip;      // Covariate
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
  real alpha4;
  real alpha5;
  real alpha6;
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
  real beta6;
  
  vector[nsites] eps;            // Random site effects
  real<lower=0> sd_eps;
  matrix[R, T] delta;            // Random transect-visit effects
  real<lower=0> sd_p;
}

transformed parameters {
  vector[R] log_lambda; // Log population size
  matrix[R, T] logit_p; // Logit detection probability

  for (i in 1:R) {
  log_lambda[i] = alpha0 + alpha1 * elev[i] + alpha2 * elev2[i] + alpha3 * twi[i] + alpha4 * litter[i] + alpha5 * gcover[i] + alpha6 * stream[i] + eps[sites[i]] * sd_eps;
  for (t in 1:T) {
  logit_p[i,t] = beta0 + beta1 * temp[i,t] + beta2 * temp2[i,t] + beta3 * precip[i,t] + beta4 * gcover[i] + beta5 * gcover2[i] + beta6 * RH[i,t] + delta[i, t] * sd_p;
  }
}
}

model {
  // Priors
  // Improper flat priors are implicitly used on alpha0, alpha1, beta0 etc. if not specified

  alpha0 ~ normal(0, 10);
  alpha1 ~ normal(0, 10);
  alpha2 ~ normal(0, 10);
  alpha3 ~ normal(0, 10);
  alpha4 ~ normal(0, 10);
  alpha5 ~ normal(0, 10);
  alpha6 ~ normal(0, 10);
  
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  beta4 ~ normal(0, 10);
  beta5 ~ normal(0, 10);
  beta6 ~ normal(0, 10);

  eps ~ normal(0, 1);
  sd_eps ~ cauchy(0, 10);
  
  for (i in 1:R) {
      for (t in 1:T) {
  delta[i,t] ~ normal(0, 1);
      }
  }
  // delta ~ normal(0, 1);
  sd_p ~ cauchy(0, 10);
  
    // Likelihood
  for (i in 1:R) {
    vector[K - max_y[i] + 1] lp;

    for (j in 1:(K - max_y[i] + 1))
      lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
             + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]);
    target += log_sum_exp(lp);
  }
  
//     for (i in 1:R) {
//     vector[K - max_y[i] + 1] lp;
// for (t in 1:T) {
//     for (j in 1:(K - max_y[i] + 1))
//       lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
//              + binomial_logit_lpmf(y[i, t] | max_y[i] + j - 1, logit_p[i, t]);
//     target += log_sum_exp(lp);
// }
//   }
}

generated quantities {
  int<lower=0> N[R];                // Abundance
  int totalN;
  vector[R] log_lik;
  // vector[N_ll] log_lik;
  // matrix[R, T] log_lik;
  // vector[K + 1] log_lik;
  // matrix[N_ll, T] log_lik;
  real mean_abundance;
  real mean_detection;
  vector[R] mean_p;
  real fit = 0;
  real fit_new = 0;
  matrix[R, T] p; 
  // matrix[R, T] p_vectorized;
  // matrix[R, T] p_test;

    matrix[R, T] eval;         // Expected values
    int y_new[R, T];
    int y_post_check[R, T];
    int y_new_sum[R];
    int y_sum_diff[R];
    matrix[R, T] y_diff;          // observed - expected
    matrix[R, T] E;
    matrix[R, T] E_new;
    //  matrix[T, 1] E[R];
    // matrix[T, 1] E_new[R];
    vector[K + 1] lp;
    // matrix[K + 1, T] lp;
   // vector[50] log_lik;
  
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
    //   for (j in 1:T) {
    //   p[i, j] = inv_logit(logit_p[i, j]);
    //   p_test[i, j] = p[i, j] - p_vectorized[i, j];
//         for (n in 0:(max_y[i] - 1))
//           lp[n + 1] = negative_infinity();
//         for (n in max_y[i]:K) {
//           lp[n + 1] = poisson_log_lpmf(n | log_lambda[i])
//             + binomial_lpmf(y[i, 1:T] | n, p[i, 1:T]);
      //   }
      // 
      mean_p[i] = mean(p[i]);
      
      for (j in 1:T) {
      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
      eval[i, j] = p[i, j] * N[i];
      E[i, j] = square(y[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
      // Generate replicate data and
      // Compute fit statistic E_new for replicate data
          y_new[i, j] = binomial_rng(N[i], p[i, j]);
          E_new[i, j] = square(y_new[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
          
          y_diff[i, j] = y[i, j] - eval[i, j];  // vectorize later observed vs. fitted
          y_post_check[i, j] = y[i, j] - y_new[i, j];
        }
        y_new_sum[i] = sum(y_new[i]);
        y_sum_diff[i] = sum(y[i]) - y_new_sum[i];
      }
    
   totalN = sum(N);  // Total pop. size across all sites
    for (i in 1:R) {
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
  mean_abundance = exp(alpha0);
  mean_detection = 1 / (1 + exp(-1 * alpha0));
}
