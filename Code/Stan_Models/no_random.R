## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.2. Generation and analysis of simulated data
## 12.2.2. Introducing covariates

library(rstan)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

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

PJORmin5 <- PJORmin[rows]
length(PJORmin5)

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
            "alpha2",
            "alpha3",
            "beta0", 
            "beta1",
            "beta2",
            "beta3",
            "beta4",
            "beta5",
            "N",
            "p",
            "log_lik",
            # "lp",
            "y_new",
            "fit",
            "fit_new")

## MCMC settings
testing <- TRUE
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


foo <- K - PJORmin5 + 1
N_ll <- sum(foo)

data.frame(num = 100:120, PJOR5[100:120, ], ymax = PJORmin5[100:120], points = foo[100:120])

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
       beta5 = runif(1, -1, 1)))

## Call Stan from R
no_random_pjor <- stan("Code/Stan_Models/no_random.stan",
            data = list(y = PJOR5, 
                        R = n.transects, 
                        T = n.surveys, 
                        elev = elev5,
                        elev2 = elev5^2,
                        litter = litter5,
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
            open_progress = FALSE, verbose = TRUE)

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
saveRDS(no_random_pjor, file = "Results/Stan/no_random_pjor_hmc.Rds")

print(no_random_pjor, digits = 3)
plot(no_random_pjor, par = c("alpha0", "alpha1", "alpha2", "beta0", "beta1"))
traceplot(no_random_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "beta4"))

print(no_random_pjor, par = "lp", digits = 2)

library("rstanarm")
library("bayesplot")
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(no_random_pjor, parameter_name = "log_lik", merge_chains = FALSE)

# removeal all -inf

# loo
r_eff <- relative_eff(exp(log_lik_1)) 
loo_no_random <- loo(log_lik_1, r_eff = r_eff, cores = nc)
print(loo_no_random)

psis_no_random <- psis(log_lik_1, r_eff = r_eff, cores = nc)
print(psis_no_random)

loo_no_random <- list(loo = loo_no_random, psis = psis_no_random, r_eff = r_eff)
saveRDS(loo_no_random, file = "Results/Stan/no_random_pjor_loo.Rds")

# Bayesian p-value check
plot(no_random_pjor, par = c("fit", "fit_new"))


