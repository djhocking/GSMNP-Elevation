############################################
# Script for checking MCMC/HMC diagnostics - this might be better as an Rmd to just have one output for supplemental publication
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
library(parallel)

#----- Set Defaults and Conditions ------

#----- Load Data and Results -----

load(file = "Data/Derived/settings.RData")
load("Data/Derived/stan_prep.RData")
site_od_full_pjor <- readRDS(file = "Results/Stan/final_od_pjor_hmc.Rds")

#----- View Traceplots (consider alt. see Simpson/Aki comments) -----

# because long enough chains almost always look "good" rank plots might be better going forward

traceplot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_eps", "sd_p"))

# mcmc_trace(site_od_full_pjor, regex_pars = "alpha")
# mcmc_trace(site_od_full_pjor, regex_pars = "N") # too many for one plot

#----- Check Domain Specific Expectations -----

# Check N for Truncation 

# The augmentation to marginalize N out as a latent discrete requires setting an upper bound, K, to loop through. If K is too small the posterior will be truncated. Need to check for every N.

# Could also set K to be so high that any truncation is either a hard prior or that the probability of detection (p) is below 10% and therefore the N-mixture model is not appopriate for those sites.

# Too many N for one plot layout so need to separate into chunks of 16 - 24

N_min <- apply(PJOR5, 1, max) # max caught on any day, so know true N is at least that large

N <- rep(NA_character_, 159)
for(i in 1:159) {
  N[i] <- paste0("N[", i, "]")
}

n <- 25

for(i in 1:(ceiling(159 / n))) {
  k <- min(i*n, 159)
  j <- k - n
  hists <- mcmc_hist(site_od_full_pjor, pars= N[j:k], binwidth = 1)
  ggsave(hists, file = paste0("Results/Stan/pjor_N_hist", i, ".pdf"))
  print(hists)
}
dev.off()
# mcmc_dens(site_od_full_pjor, pars = c("N[100]", "N[120]"))

# check detection
# In simulations abundance esimates are unreliable unless detection > 0.20

mcmc_dens(site_od_full_pjor, pars = c("mean_detection"))

#----- Check Divergences and Pairwise Correlations -----

check_divergences(site_od_full_pjor)

color_scheme_set("darkgray")
# mcmc_scatter(
#   as.matrix(site_od_full_pjor),
#   pars = c("fit", "fit_new"), 
#   np = nuts_params(site_od_full_pjor), 
#   np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
# )

# mcmc_scatter(
#   as.matrix(site_od_full_pjor),
#   pars = c("eval[1,1]", "y_new[1,1]"), 
#   np = nuts_params(site_od_full_pjor), 
#   np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
# )

# look at variance terms because usually the hardest to fit
color_scheme_set("darkgray")
# mcmc_scatter(
#   site_od_full_pjor,
#   pars = c("sd_eps", "sd_p"), 
#   np = nuts_params(site_od_full_pjor), 
#   np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
# )

mcmc_pairs(
  site_od_full_pjor,
  pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p", "lp__"), 
  off_diag_args = list(alpha = 0.1),
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

# pairs(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p", "lp__"))

#----- Check Energy and Treedepth -----

check_energy(site_od_full_pjor)

check_treedepth(site_od_full_pjor)


#----- Summarize Samples Sizes and Mixing -----

# Effective samples sizes using rstan::monitor following Hoffman and Gelman (2014) to be more reliable and accurate (as in Monnahan et al. 2017) - UPDATE - now follow Vehtari et al. 2019.

# effective sample sizes (should be > 100) - Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2019). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
mon <- rstan::monitor(site_od_full_pjor, warmup = nb)
# rstan:::print.simsummary(mon)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

#----- Posterior Predictive Checks -----


# Bayesian p-value check
plot(site_od_full_pjor, par = c("fit", "fit_new"))
pairs(site_od_full_pjor, pars = c("fit", "fit_new"))

fit <- rstan::extract(site_od_full_pjor, par = "fit")[[1]]
fit_new <- rstan::extract(site_od_full_pjor, par = "fit_new")[[1]]
# plot(fit, fit_new)

mean(fit_new > fit) # Bayesian p-value - not sure why it's much worse than JAGS output - especially in light of posterior predictive checks being great below


# examine posterior prredictions of total counts across all 5 visits
y_sum <- data.frame(transect = rownames(PJOR5), y_sum = rowSums(PJOR5), stringsAsFactors = FALSE)
y_sum_new <- rstan::extract(site_od_full_pjor, pars = c("y_new_sum"))

ppc_scatter_avg(y = y_sum$y_sum, yrep = y_sum_new[[1]]) # average posterior across samples

ppc_scatter(y = y_sum$y_sum, yrep = y_sum_new[[1]][sample(1:nrow(y_sum_new[[1]]), size = 16, replace = FALSE) , ]) # iteration specific for 16 random otherwise too hard to plot. Maybe could overlay but probably unnecessary to see all iterations when things generally look very good

# RMSE of posterior predictive
sqrt(sum(((y_sum$y_sum - colMeans(y_sum_new[[1]]))^2) / nrow(y_sum)))

# Posterior predictive check for each visit
y_long <- PJOR5 %>%
  mutate(site = rownames(.)) %>%
  pivot_longer(starts_with("PJOR"),
               names_to = "visit",
               names_prefix = "PJOR")

y_new<- rstan::extract(site_od_full_pjor, pars = c("y_new"))[[1]]

y_new_mat <- matrix(NA_integer_, dim(y_new)[1], nrow(y_long))
for(i in 1:dim(y_new)[1]) {
  y_new_long <- y_new[i, , ] %>%
    data.frame() %>%
    mutate(site = rownames(.)) %>%
    pivot_longer(starts_with("X"),
                 names_to = "visit",
                 names_prefix = "X")
  y_new_mat[i, ] <- y_new_long$value
}

# RMSE of posterior predictive for observations per visit
sqrt(sum(((y_long$value - colMeans(y_new_mat))^2) / nrow(y_sum)))

ppc_scatter_avg(y = y_long$value, yrep = y_new_mat, alpha = 0.3)

ppc_rootogram(y = y_long$value, yrep = y_new_mat, style = "standing")
# ppc_rootogram(y = y_long$value, yrep = y_new_mat, style = "hanging")
# ppc_rootogram(y = y_long$value, yrep = y_new_mat, stayle = "suspended")
ppc_bars(y = y_long$value, yrep = y_new_mat)
ppc_bars_grouped(y = y_long$value, yrep = y_new_mat, group = y_long$visit)


# RMSE of posterior predictive
sqrt(sum(((y_long$value - colMeans(y_new_mat))^2) / nrow(y_long)))


#----- Prior vs. Postior Distributions -----


#----- Leave one (site) out cross validation -----

# print(site_od_full_pjor, par = "log_lik", digits = 2)

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(site_od_full_pjor, parameter_name = "log_lik", merge_chains = FALSE)

# remove all -inf

# loo
r_eff <- relative_eff(exp(log_lik_1)) 
loo_od <- loo::loo(log_lik_1, r_eff = r_eff, cores = nc)
print(loo_od)
pareto_k_table(loo_od)

loo_od <- loo::loo(log_lik_1, r_eff = r_eff, cores = parallel::detectCores(), K = 10)
loo_od
plot(loo_od)

# yrep <- posterior_predict(site_od_full_pjor)
# 
# ppc_loo_pit_overlay(
#   y = roaches$y,
#   yrep = yrep,
#   lw = weights(loo1$psis_object)
# )

psis_od <- psis(log_lik_1, r_eff = r_eff, cores = nc)
print(psis_od)
pareto_k_table(psis_od)
plot(psis_od, label_points = TRUE)

loo_od <- list(loo = loo_od, psis = psis_od, r_eff = r_eff)
saveRDS(loo_od, file = "Results/Stan/site_od_full_pjor_loo.Rds")


#----- Other -----

list_of_draws <- rstan::extract(site_od_full_pjor)
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
