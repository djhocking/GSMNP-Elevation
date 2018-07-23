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
PJOR5_na <- PJOR[ , 1:5]

na_rows <- which(is.na(rowSums(PJOR5_na)))
rows <- which(!is.na(rowSums(PJOR5_na)))
PJOR5 <- PJOR5_na[which(!is.na(rowSums(PJOR5_na))), ]
dim(PJOR5)
summary(PJOR5)

elev5 <- elev[rows]
stream5 <- stream[rows]
litter5 <- litterdepth[rows]

RH5 <- RH.s[rows, 1:5]
temp5 <- Temp.s[rows, 1:5]
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
            # "alpha2",
            "alpha3",
            "beta0", 
            "beta1",
            "beta2",
            "beta3",
            "beta4",
            # "beta5",
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
  nb = 400
  ni = 600
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
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
       # beta5 = runif(1, -1, 1),
       sd_eps = runif(1, 0, 1)))

## Call Stan from R
if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
site_od_elev_pjor <- stan("Code/Stan_Models/site_od_elev.stan",
                     data = list(y = PJOR5, 
                                 R = n.transects, 
                                 T = n.surveys, 
                                 nsites = n.sites,
                                 sites = Data5$site_stan,
                                 elev = elev5,
                                 # elev2 = elev5^2,
                                 litter = litter5,
                                 gcover = gcover5,
                                 # gcover2 = gcover5^2,
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
saveRDS(site_od_elev_pjor, file = "Results/Stan/site_od_elev_pjor_hmc.Rds")

print(site_od_elev_pjor, digits = 3)
plot(site_od_elev_pjor, par = c("alpha0", "alpha1", "beta0", "beta1", "fit", "fit_new"))
traceplot(site_od_elev_pjor, par = c("alpha0", "alpha1", "alpha3", "beta0", "beta1", "beta2", "beta3", "beta4", "sd_eps", "sd_p"))

pairs(site_od_elev_pjor, pars = c("alpha0", "alpha1", "beta0", "sd_eps", "sd_p"))

# print(site_od_elev_pjor, par = "log_lik", digits = 2)

library("rstanarm")
library("bayesplot")
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_2 <- extract_log_lik(site_od_elev_pjor, parameter_name = "log_lik", merge_chains = FALSE)

# removeal all -inf

# as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- relative_eff(exp(log_lik_2)) 
loo_2 <- loo(log_lik_2, r_eff = r_eff, cores = nc)
print(loo_2)

psis_od_elev <- psis(log_lik_2)

# Bayesian p-value check
plot(site_od_elev_pjor, par = c("fit", "fit_new"))

foo <- as.data.frame(summary(site_od_elev_pjor)$summary)
foo$parameter <- rownames(summary(site_od_elev_pjor)$summary)

loo::compare(loo_1, loo_2)

# waic1 <- waic(log_lik_1)
# waic2 <- waic(log_lik_2)
# compare(waic1, waic2)

library(stringr)
y_new_sum <- dplyr::filter(foo, str_detect(parameter, "y_new"))
eval_sum <- dplyr::filter(foo, str_detect(parameter, "eval"))

plot(PJOR5, y_new_sum$mean)
abline(a = 0, b = 1)

plot(eval_sum, y_new_sum$mean)
abline(a = 0, b = 1)
