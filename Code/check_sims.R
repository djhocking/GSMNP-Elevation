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

site_od_full_pjor <- readRDS(file = "Results/Stan/final_od_pjor_hmc.Rds")

#----- View Traceplots (consider alt. see Simpson/Aki comments) -----

# because long enough chains almost always look "good" rank plots might be better going forward

traceplot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_eps", "sd_p"))

# mcmc_trace(site_od_full_pjor, regex_pars = "alpha")
# mcmc_trace(site_od_full_pjor, regex_pars = "N") # too many for one plot

#----- Check Divergences and Pairwise Correlations -----

check_divergences(site_od_full_pjor)

color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(site_od_full_pjor),
  pars = c("y_diff[1,1]", "y_diff[1,2]"), 
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

mcmc_scatter(
  as.matrix(site_od_full_pjor),
  pars = c("fit", "fit_new"), 
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

mcmc_scatter(
  as.matrix(site_od_full_pjor),
  pars = c("eval[1,1]", "y_new[1,1]"), 
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(site_od_full_pjor),
  pars = c("sd_eps", "sd_p"), 
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

mcmc_pairs(
  as.matrix(site_od_full_pjor),
  pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p"), 
  np = nuts_params(site_od_full_pjor), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

pairs(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "beta0", "sd_eps", "sd_p", "lp__"))

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
plot(fit, fit_new)

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
  pivot_longer(starts_with("P"),
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

ppc_scatter_avg(y = y_long$value, yrep = y_new_mat, alpha = 0.3)

ppc_rootogram(y = y_long$value, yrep = y_new_mat, style = "standing")
# ppc_rootogram(y = y_long$value, yrep = y_new_mat, style = "hanging")
# ppc_rootogram(y = y_long$value, yrep = y_new_mat, stayle = "suspended")
ppc_bars(y = y_long$value, yrep = y_new_mat)
ppc_bars_grouped(y = y_long$value, yrep = y_new_mat, group = y_long$visit)


# RMSE of posterior predictive
sqrt(sum(((y_long$value - colMeans(y_new_mat))^2) / nrow(y_long)))


#----- Prior vs. Postior Distributions -----



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
