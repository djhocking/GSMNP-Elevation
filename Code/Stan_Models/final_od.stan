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
  int<lower=0> K[R];       // Upper bound of population size
}

transformed data {
  int<lower=0> max_y[R];
  int<lower=0> N_ll;
  int tmp[R];
  
  for (i in 1:R) {
    max_y[i] = max(y[i]);
    tmp[i] = K[i] - max_y[i] + 1;
  }
  N_ll = sum(tmp);
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
    log_lambda[i] = alpha0 + alpha1 * elev[i] + alpha2 * elev2[i] + alpha3 * twi[i] + alpha4 * litter[i] + alpha5 * gcover[i] + alpha6 * stream[i] + eps[sites[i]] * sd_eps; // non-centered formulation of random effect (see Monnahan et al. 2017)
    for (t in 1:T) {
      logit_p[i,t] = beta0 + beta1 * temp[i,t] + beta2 * temp2[i,t] + beta3 * precip[i,t] + beta4 * gcover[i] + beta5 * gcover2[i] + beta6 * RH[i,t] + delta[i, t] * sd_p; // non-centered formulation
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
  
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 2);
  beta2 ~ normal(0, 2);
  beta3 ~ normal(0, 2);
  beta4 ~ normal(0, 2);
  beta5 ~ normal(0, 2);
  beta6 ~ normal(0, 2);
  
  eps ~ normal(0, 1);
  sd_eps ~ cauchy(0, 2.5);
  
  for (i in 1:R) {
    for (t in 1:T) {
      delta[i,t] ~ normal(0, 1);
    }
  }
  
  sd_p ~ cauchy(0, 2.5); // even this might be more heavily tailed than we want on p. maybe half normal with sd = 3-5?
  
  // Likelihood
  for (i in 1:R) {
    vector[K[i] - max_y[i] + 1] lp;
    
    for (j in 1:(K[i] - max_y[i] + 1))
    lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
    + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]); // vectorized over T
    target += log_sum_exp(lp);
  }
}

generated quantities {
  // int<lower=0> N[R];                  // Abundance (must be at least max_y)
  int N[R] = max_y;
  // int<lower=0> N_rng[R];              // Abundance proposed by random draw from lambda
  int N_total;
  vector[R] log_lik;
  real mean_abundance;
  real mean_detection;
  real mean_p;
  vector[R] mean_p_site;
  real fit = 0;
  real fit_new = 0;
  matrix[R, T] p; 
  matrix[R, T] eval;         // Expected values
  int y_new[R, T];
  int y_new_sum[R];
  matrix[R, T] E;
  matrix[R, T] E_new;
  int counter[R];
  // real<lower=0> lambda[R];
  
  // for(i in 1:R) {
  //   
  // }
  
  for (i in 1:R) {
    // real sum_p = 0;
    // real u = uniform_rng(0, 1);
    
    // calculate vector logliklihood for use in loo for model comparison
    vector[K[i] - max_y[i] + 1] ll;
    
    for (k in 1:(K[i] - max_y[i] + 1)) {
      ll[k] = poisson_log_lpmf(max_y[i] + k - 1 | log_lambda[i]) // remake lp because can't use from model statement in generated quantities
      + binomial_logit_lpmf(y[i] | max_y[i] + k - 1, logit_p[i]);
    }
    log_lik[i] = log_sum_exp(ll); // for use in loo and multimodel comparison

    // Calculate Abundance - Restrict N to be at least as big as the number of animals observed on a single night
      N[i] = poisson_log_rng(log_lambda[i]);
      counter[i] = 0;
      while (N[i] < max_y[i]) {
        N[i] = poisson_log_rng(log_lambda[i]);
        counter[i] += 1;
        if (counter[i] > 100) break;
      }
      
      // N[i] = poisson_log_rng(log_lambda[i]);
      // lambda[i] = exp(log_lambda[i]);
      // 
      // for (k in max_y[i]:K[i]) {
      //   sum_p = sum_p + exp(poisson_lpmf(k | lambda[i]) - poisson_lcdf(K[i] | lambda[i]));
      //   if (sum_p >= u) {
      //     N[i] = k;
      //     break;
      //   }
      // }
      
    p[i, 1:T] = inv_logit(logit_p[i, 1:T]);
  }
  
  // Bayesian p-value fit - Diff than JAGS but all parameters about the same. Not a great test anyway. Problem with RNG for poisson rather than sampling N directly. Makes a big diff with the chi-sq and other Bayesian P-value metrics. But other posterior predictive checks look great.
  
  // Initialize E and E_new
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
    mean_p_site[i] = mean(p[i]);
    
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
    y_new_sum[i] = sum(y_new[i]);
    fit = fit + sum(E[i]);
    fit_new = fit_new + sum(E_new[i]);
  }
  
  N_total = sum(N);  // Total pop. size across all sites
  mean_abundance = exp(alpha0);
  mean_detection = 1 / (1 + exp(-1 * beta0));
  mean_p = mean(p);
}
