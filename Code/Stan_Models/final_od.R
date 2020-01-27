## Random site effects on abundance and random binomial overdispersion with random transect-visit effect

# Based on code by Hiroki ITÔ translated from chapter 12 of Kery and Schuab Bayesian Population Analysis 

# non-center as described by Monnahan et al. 2017 MEE

# Load Libraries
library(rstan)
library(dplyr)

# Settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set.seed(123)

load(file = "Data/Derived/settings.RData")

## Read data
load("Data/Derived/stan_prep.RData")


## Parameters monitored
params <- c("N_total", 
            "alpha0", 
            "alpha1", 
            "alpha2",
            "alpha3",
            "alpha4",
            "alpha5",
            "alpha6",
            "beta0", 
            "beta1",
            "beta2",
            "beta3",
            "beta4",
            "beta5",
            "beta6",
            "sd_eps", 
            "sd_p",
            "N",
            # "p",
            "mean_abundance",
            "mean_detection",
            "mean_p",
            "log_lik",
            # "p_test",
            # "lp",
            "eval",
            "y_new",
            # "y_diff",
            # "y_post_check",
            "y_new_sum",
            # "y_sum_diff",
            "fit",
            "fit_new")



## Initial values
inits <- lapply(1:nc, function(i)
  list(alpha0 = runif(1, -1, 1),
       alpha1 = runif(1, -1, 1),
       alpha2 = runif(1, -1, 1),
       alpha3 = runif(1, -1, 1),
       alpha4 = runif(1, -1, 1),
       alpha5 = runif(1, -1, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
       beta5 = runif(1, -1, 1),
       sd_eps = runif(1, 0, 1)))

## Call Stan from R

#----- EWIL -----
if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
site_od_full_ewil <- stan("Code/Stan_Models/final_od.stan",
                          data = list(y = EWIL5, 
                                      R = nrow(EWIL5), 
                                      T = ncol(EWIL5), 
                                      nsites = n.sites,
                                      sites = Data5$site_stan,
                                      elev = elev5,
                                      elev2 = elev5^2,
                                      litter = litter5,
                                      twi = twi5,
                                      precip = precip5,
                                      stream = stream5,
                                      # stream2 = stream5 * stream5,
                                      gcover = gcover5,
                                      gcover2 = gcover5^2,
                                      RH = RH5,
                                      temp = temp5,
                                      temp2 = temp5^2,
                                      K = K_ewil),
                          init = inits,
                          pars = params,
                          chains = nc, iter = ni, warmup = nb, thin = nt,
                          # seed = 1,
                          open_progress = FALSE, 
                          verbose = TRUE)

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
saveRDS(site_od_full_ewil, file = "Results/Stan/final_od_ewil_hmc.Rds")


#----- DWRI -----

site_od_full_dwri <- stan("Code/Stan_Models/final_od.stan",
                          data = list(y = DWRI5, 
                                      R = nrow(DWRI5), 
                                      T = ncol(DWRI5), 
                                      nsites = n.sites,
                                      sites = Data5$site_stan,
                                      elev = elev5,
                                      elev2 = elev5^2,
                                      litter = litter5,
                                      twi = twi5,
                                      precip = precip5,
                                      stream = stream5,
                                      # stream2 = stream5 * stream5,
                                      gcover = gcover5,
                                      gcover2 = gcover5^2,
                                      RH = RH5,
                                      temp = temp5,
                                      temp2 = temp5^2,
                                      K = K_dwri),
                          init = inits,
                          pars = params,
                          chains = nc, iter = ni, warmup = nb, thin = nt,
                          # seed = 1,
                          open_progress = FALSE, 
                          verbose = TRUE)

saveRDS(site_od_full_dwri, file = "Results/Stan/final_od_dwri_hmc.Rds")


#----- PJOR -----

site_od_full_pjor <- stan("Code/Stan_Models/final_od.stan",
                          data = list(y = PJOR5, 
                                      R = nrow(PJOR5), 
                                      T = ncol(PJOR5), 
                                      nsites = n.sites,
                                      sites = Data5$site_stan,
                                      elev = elev5,
                                      elev2 = elev5^2,
                                      litter = litter5,
                                      twi = twi5,
                                      precip = precip5,
                                      stream = stream5,
                                      # stream2 = stream5 * stream5,
                                      gcover = gcover5,
                                      gcover2 = gcover5^2,
                                      RH = RH5,
                                      temp = temp5,
                                      temp2 = temp5^2,
                                      K = K_pjor),
                          init = inits,
                          pars = params,
                          chains = nc, iter = ni, warmup = nb, thin = nt,
                          # seed = 1,
                          open_progress = FALSE, 
                          verbose = TRUE)
saveRDS(site_od_full_pjor, file = "Results/Stan/final_od_pjor_hmc.Rds")

#---- Leftovers -----

# print(site_od_full_pjor, par = "log_lik", digits = 2)

library("rstanarm")
library("bayesplot")
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(site_od_full_pjor, parameter_name = "log_lik", merge_chains = FALSE)

# remove all -inf

# loo
r_eff <- relative_eff(exp(log_lik_1)) 
loo_od <- loo::loo(log_lik_1, r_eff = r_eff, cores = nc)
print(loo_od)

psis_od <- psis(log_lik_1, r_eff = r_eff, cores = nc)
print(psis_od)

loo_od <- list(loo = loo_od, psis = psis_od, r_eff = r_eff)
saveRDS(loo_od, file = "Results/Stan/site_od_full_pjor_loo.Rds")

#----- Cleanup -----

rm(list = ls())

