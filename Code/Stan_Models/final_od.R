## Random site effects on abundance and random binomial overdispersion with random transect-visit effect

# Based on code by Hiroki ITÔ translated from chapter 12 of Kery and Schuab Bayesian Population Analysis 

# non-center as described by Monnahan et al. 2017 MEE

# Load Libraries
library(rstan)
library(dplyr)

# Settings
testing <- FALSE # settings to run quickly when testing model changes = TRUE

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set.seed(123)

## Read data
load("Data/Processed/jags_prep.RData")

## Remove missing observations for now
summary(PJOR)
PJOR5_na <- PJOR[ , 1:5]

na_rows <- which(is.na(rowSums(PJOR5_na)))
rows <- which(!is.na(rowSums(PJOR5_na)))
PJOR5 <- PJOR5_na[which(!is.na(rowSums(PJOR5_na))), ]
dim(PJOR5)
summary(PJOR5)

elev5 <- elev[rows]
stream5 <- stream[rows]
twi5 <- twi[rows]
litter5 <- litterdepth[rows]

RH5 <- RH.s[rows, 1:5]
temp5 <- Temp.s[rows, 1:5]
precip5 <- Precip.s[rows, 1:5]
gcover5 <- gcover[rows]

n.transects <- length(elev5)
n.surveys <- ncol(PJOR5)

Data5 <- Data[rows, ]
Data5 <- Data5 %>%
  arrange(site, transect)
Data5$site_stan <- 99999
Data5$site_stan[1] <- 1
for(i in 2:n.transects) {
  Data5$site_stan[i] <- ifelse(Data5$site[i] == Data5$site[i-1], Data5$site_stan[i-1], Data5$site_stan[i-1] + 1)
}

n.sites <- length(unique(Data5$site_stan))

## Parameters monitored
params <- c("totalN", 
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
            "p",
            "mean_abundance",
            "mean_detection",
            "mean_p",
            "log_lik",
            # "p_test",
            # "lp",
            "eval",
            "y_new",
            "y_diff",
            "y_post_check",
            "y_new_sum",
            "y_sum_diff",
            "fit",
            "fit_new")

## MCMC settings
if(testing) {
  nb = 100
  ni = 200
  nt = 1
  nc <- 3
  K = max(apply(PJOR5, MARGIN = 1, max)) + 1
} else {
  nb = 1000
  ni = 1500
  nt = 4
  nc <- parallel::detectCores()
  K = 500 # should go back to getting variable K code working even if more bookeeping - waiting massive time and resources on 90% of sites
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
                                      K = K),
                          init = inits,
                          pars = params,
                          chains = nc, iter = ni, warmup = nb, thin = nt,
                          # seed = 1,
                          open_progress = FALSE, 
                          verbose = TRUE)

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
saveRDS(site_od_full_pjor, file = "Results/Stan/final_od_pjor_hmc.Rds")

# effective sample sizes (should be > 100) - Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2019). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
mon <- monitor(site_od_full_pjor, warmup = nb)
rstan:::print.simsummary(mon)
# print(mon)
# print(site_od_full_pjor)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

# print(site_od_full_pjor, digits = 3)
# print(site_od_full_pjor, pars = "p_test", digits = 3)
plot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6"))
plot(site_od_full_pjor, par = c("y_diff"))
traceplot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_eps", "sd_p"))

library(bayesplot)
mcmc_trace(site_od_full_pjor, regex_pars = "N")

plot(site_od_full_pjor, par = "y_sum_diff")

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


# divergences
library(bayesplot)

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

mean(fit_new > fit) # Bayesian p-value - not sure why it's much worse than JAGS output - maybe because of summing with NA being handled differently

foo <- as.data.frame(summary(site_od_full_pjor)$summary)
foo$parameter <- rownames(summary(site_od_full_pjor)$summary)

library(stringr)
y_new_sum <- dplyr::filter(foo, str_detect(parameter, "y_new"))
eval_sum <- dplyr::filter(foo, str_detect(parameter, "eval"))

plot(as.numeric(unlist(PJOR5)), y_new_sum$mean)
abline(a = 0, b = 1)

plot(eval_sum$mean, y_new_sum$mean)
abline(a = 0, b = 1)


#----- Posterior Predictive Checks -----
library(tidyr)

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
