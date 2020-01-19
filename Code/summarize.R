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

#----- Summarize coefficients -----

#----- Plot parameter summaries -----

plot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6"))
plot(site_od_full_pjor, par = c("y_diff"))

# plot(site_od_full_pjor, par = "y_sum_diff")
