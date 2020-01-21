
load(file = "Data/Derived/stan_prep.RData")

if(!dir.exists("Data/Derived")) dir.create("Data/Derived", recursive = TRUE)

testing <- FALSE # settings to run quickly when testing model changes = TRUE

## MCMC settings
if(testing) {
  nb = 100
  ni = 200
  nt = 1
  nc <- 3
  K_pjor = apply(PJOR5, 1, max) + 10
  K_ewil = apply(EWIL5, 1, max) + 10
  K_dwri = apply(DWRI5, 1, max) + 10
} else {
  nb = 1000
  ni = 2000
  nt = 2
  nc <- 6 # parallel::detectCores()
  K_pjor = as.integer((apply(PJOR5, 1, max)) / 0.1 + 10) # should go back to getting variable K code working even if more bookeeping - waiting massive time and resources on 90% of sites
  K_ewil = as.integer((apply(EWIL5, 1, max)) / 0.1 + 10)
  K_dwri = as.integer((apply(DWRI5, 1, max)) / 0.1 + 10)
}

save(testing, nb, ni, nt, nc, K_pjor, K_ewil, K_dwri, file = "Data/Derived/settings.RData")

rm(list = ls())