## Random site effects on abundance and random binomial overdispersion with random transect-visit effect

# Based on code by Hiroki ITÔ translated from chapter 12 of Kery and Schuab Bayesian Population Analysis 

# Load Libraries
library(rstan)
library(dplyr)

# Settings
testing <- TRUE # settings to run quickly when testing model changes = TRUE

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set.seed(123)

## Read data
load("Data/Processed/jags_prep.RData")

## Remove missing observations for now
summary(PJOR)

visits <- rowSums(!is.na(PJOR))

## Parameters monitored
params <- c("totalN", 
            "alpha0", 
            "alpha1", 
            "alpha2",
            "alpha3",
            "alpha4",
            "alpha5",
            "beta0", 
            "beta1",
            "beta2",
            "beta3",
            "beta4",
            "beta5",
            "sd_eps", 
            "sd_p",
            "N",
            "p",
            "log_lik",
            # "lp",
            "eval",
            "y_new",
            "fit",
            "fit_new")

## MCMC settings
if(testing) {
  nb = 300
  ni = 400
  nt = 1
  nc <- 3
  K = 100
} else {
  nb = 5000
  ni = 10000
  nt = 1
  nc <- 4
  K = 300
}

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
if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
site_od_full_pjor <- stan("Code/Stan_Models/site_od_full.stan",
                          data = list(y = PJOR5, 
                                      R = n.transects, 
                                      T = n.surveys, 
                                      nsites = n.sites,
                                      sites = Data5$site_stan,
                                      elev = elev5,
                                      elev2 = elev5^2,
                                      litter = litter5,
                                      stream = stream5,
                                      stream2 = stream5 * stream5,
                                      gcover = gcover5,
                                      gcover2 = gcover5^2,
                                      RH = RH5,
                                      temp = temp5,
                                      temp2 = temp5^2,
                                      K = K),
                          init = inits,
                          pars = params,
                          chains = nc, iter = ni, warmup = nb, thin = nt,
                          # seed = 1,
                          open_progress = FALSE, 
                          verbose = TRUE)

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
saveRDS(site_od_full_pjor, file = "Results/Stan/site_od_full_pjor_hmc.Rds")

print(site_od_full_pjor, digits = 3)
plot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5"))
traceplot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "beta4", "sd_eps", "sd_p"))

pairs(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p", "lp__"))

list_of_draws <- extract(site_od_full_pjor)
print(names(list_of_draws))

print(get_elapsed_time(site_od_full_pjor))

sampler_params <- get_sampler_params(site_od_full_pjor, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

str(sampler_params)

sapply(sampler_params, function(x) mean(x[, "energy__"]))

energy__1 = as.data.frame(sampler_params[[1]])$energy__
alpha0 <- unlist(site_od_full_pjor@sim$samples[[1]]["alpha0"])[(nb+1):ni]
alpha1 <- unlist(site_od_full_pjor@sim$samples[[1]]["alpha1"])[(nb+1):ni]
alpha2 <- unlist(site_od_full_pjor@sim$samples[[1]]["alpha2"])[(nb+1):ni]
beta0 <- unlist(site_od_full_pjor@sim$samples[[1]]["beta0"])[(nb+1):ni]
sd_eps <- unlist(site_od_full_pjor@sim$samples[[1]]["sd_eps"])[(nb+1):ni]
sd_p <- unlist(site_od_full_pjor@sim$samples[[1]]["sd_p"])[(nb+1):ni]
lp <- unlist(site_od_full_pjor@sim$samples[[1]]["lp__"])[(nb+1):ni]
pairs(data.frame(energy__1, alpha0, alpha1, alpha2, beta0, sd_eps, sd_p))


mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)


print(site_od_full_pjor, par = "log_lik", digits = 2)

library("rstanarm")
library("bayesplot")
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(site_od_full_pjor, parameter_name = "log_lik", merge_chains = FALSE)

# remove all -inf

# loo
r_eff <- relative_eff(exp(log_lik_1)) 
loo_od <- loo(log_lik_1, r_eff = r_eff, cores = nc)
print(loo_od)

psis_od <- psis(log_lik_1, r_eff = r_eff, cores = nc)
print(psis_od)

loo_od <- list(loo = loo_od, psis = psis_od, r_eff = r_eff)
saveRDS(loo_od, file = "Results/Stan/site_od_full_pjor_loo.Rds")

# Bayesian p-value check
plot(site_od_full_pjor, par = c("fit", "fit_new"))

pairs(site_od_full_pjor, pars = c("fit", "fit_new"))

fit <- c(site_od_full_pjor@sim$samples[[1]]$fit[(nb+1):ni], site_od_full_pjor@sim$samples[[2]]$fit[(nb+1):ni], site_od_full_pjor@sim$samples[[3]]$fit[(nb+1):ni]) # maybe do with lapply?

fit_new <- c(site_od_full_pjor@sim$samples[[1]]$fit_new[(nb+1):ni], site_od_full_pjor@sim$samples[[2]]$fit_new[(nb+1):ni], site_od_full_pjor@sim$samples[[3]]$fit_new[(nb+1):ni]) # maybe do with lapply?

mean(fit_new > fit) # Bayesian p-value - not sure why it's much worse than JAGS output

foo <- as.data.frame(summary(site_od_full_pjor)$summary)
foo$parameter <- rownames(summary(site_od_full_pjor)$summary)

library(stringr)
y_new_sum <- dplyr::filter(foo, str_detect(parameter, "y_new"))
eval_sum <- dplyr::filter(foo, str_detect(parameter, "eval"))

plot(PJOR5, y_new_sum$mean)
abline(a = 0, b = 1)

plot(eval_sum, y_new_sum$mean)
abline(a = 0, b = 1)
