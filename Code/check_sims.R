############################################
# Script for checking MCMC/HMC diagnostics
############################################

#----- Load Libraries ------

library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
# library(tidybayes)
library(loo)
# library(rstanarm)
library(bayesplot)

#----- Set Defaults and Conditions ------

#----- Load Data and Results -----


#----- View Traceplots (consider alt. see Simpson/Aki comments) -----

# because long enough chains almost always look "good" rank plots might be better going forward

traceplot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_eps", "sd_p"))

# mcmc_trace(site_od_full_pjor, regex_pars = "alpha")

# mcmc_trace(site_od_full_pjor, regex_pars = "N") # too many for one plot

#----- Check Divergences -----

#----- Check Energy Correlations -----

#----- Visualize Pairwise Correlations -----

# intercepts (alpha0 and beta0) and random standard deviations (sd_*) tend to be the most difficult to separate with poor mixing, high autocorrelation, and high correlation with the others

pairs(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p", "lp__"))

#----- Summarize Samples Sizes and Mixing -----

# effective sample sizes (should be > 100) - Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2019). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
mon <- monitor(site_od_full_pjor, warmup = nb)
# rstan:::print.simsummary(mon)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

#----- Posterior Predictive Checks -----


#----- Prior vs. Postior Distributions -----


