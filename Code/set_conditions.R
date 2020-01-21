
load(file = "Data/Derived/stan_prep.RData")

if(!dir.exists("Data/Derived")) dir.create("Data/Derived", recursive = TRUE)

testing <- TRUE # settings to run quickly when testing model changes = TRUE

## MCMC settings
if(testing) {
  nb = 100
  ni = 200
  nt = 1
  nc <- 3
  K = apply(PJOR5, 1, max) + 10
} else {
  nb = 1000
  ni = 1500
  nt = 4
  nc <- parallel::detectCores()
  K = as.integer((apply(PJOR5, 1, max)) / 0.1 + 10) # should go back to getting variable K code working even if more bookeeping - waiting massive time and resources on 90% of sites
}

save(testing, nb, ni, nt, nc, K, file = "Data/Derived/settings.RData")

rm(list = ls())