## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.2. Generation and analysis of simulated data
## 12.2.2. Introducing covariates

library(rstan)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# set.seed(123)

## Read data
load("Data/Processed/jags_prep.RData")

## Remove missing observations for now
summary(DWRI)
DWRI5_na <- DWRI[ , 1:5]

na_rows <- which(is.na(rowSums(DWRI5_na)))
rows <- which(!is.na(rowSums(DWRI5_na)))
DWRI5 <- DWRI5_na[which(!is.na(rowSums(DWRI5_na))), ]
dim(DWRI5)
summary(DWRI5)

elev5 <- elev[rows]
stream5 <- stream[rows]
litter5 <- litterdepth[rows]

RH5 <- RH.s[rows, 1:5]
temp5 <- Temp.s[rows, 1:5]
gcover5 <- gcover[rows]

n.transects <- length(elev5)
n.surveys <- ncol(DWRI5)

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
            "sd_eps", 
            "sd_p", 
            "N")

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
       beta5 = runif(1, -1, 1),
       sd_eps = runif(1, 0, 1)))

## Call Stan from R
out <- stan("Code/Stan_Models/site_od.stan",
            data = list(y = DWRI5, 
                        R = n.transects, 
                        T = n.surveys, 
                        elev = elev5,
                        litter = litter5,
                        gcover = gcover5,
                        RH = RH5,
                        temp = temp5,
                        nsites = n.sites,
                        sites = Data5$site_stan,
                        K = K),
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE, verbose = TRUE)

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
saveRDS(out, file = "Results/Stan/site_od_dwri_hmc.Rds")

print(out, digits = 3)
plot(out, par = c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "sd_eps"))
traceplot(out, par = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "sd_eps", "sd_p"))

library("rstanarm")
library("bayesplot")
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(fit_1, merge_chains = FALSE)

# as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- relative_eff(exp(log_lik_1)) 

loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
