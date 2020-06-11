# Compare models with loo
# Leave-one-site-out loo-cv

# load libraries
library(rstan)
library(loo)

# PJOR
# Load loo results
loo_od_elev_pjor <- readRDS("Results/Stan/site_od_elev_pjor_loo.Rds")
loo_od_pjor <- readRDS("Results/Stan/site_od_pjor_loo.Rds")
loo_od_full_pjor <- readRDS("Results/Stan/site_od_full_pjor_loo.Rds")
loo_od_original_pjor <- readRDS("Results/Stan/site_od_original_pjor_loo.Rds")
loo_no_random_pjor <- readRDS("Results/Stan/no_random_pjor_loo.Rds")

# compare models
loo::loo_compare(list(
  no_random = loo_no_random_pjor$loo,
  elev2_simple = loo_od_pjor$loo,
  elev_simple = loo_od_elev_pjor$loo,
  elev2 = loo_od_full_pjor$loo,
  original = loo_od_original_pjor$loo))
