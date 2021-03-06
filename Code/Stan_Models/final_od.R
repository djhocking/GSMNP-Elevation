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
            # "mean_p_site",
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


#----- Cleanup -----

rm(list = ls())

