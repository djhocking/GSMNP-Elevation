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
if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)
out <- stan("Code/Stan_Models/site_od.stan",
            data = list(y = PJOR5, 
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

print(out, digits = 3)
plot(out, par = c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "sd_eps"))
traceplot(out, par = c("alpha0", "alpha1", "alpha2", "alpha3", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "sd_eps", "sd_p"))
