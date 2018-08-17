# Summarize JAGS results

# Load libraries
library(rjags)
library(coda)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load Data
load("Data/Processed/jags_prep.RData")
load(file = "Results/JAGS/pjor_mcmc_out.RData")
load(file = "Results/JAGS/dwri_mcmc_out.RData")
load(file = "Results/JAGS/ewil_mcmc_out.RData")


##### PJOR #####
# Check fit
for(i in 1:3) bayesP.pjor4 <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.353 good
print(bayesP.pjor4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), dig=3) # good convergence

Quants.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

pjor.variables <- c("N-intercept", "Elevation", "TPI", "TWI", "Canopy", "Litter Depth", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

pjor.summary <- data.frame(pjor.variables, Means.pjor, SDs.pjor, Quants.pjor["2.5%", ], Quants.pjor["50%", ], Quants.pjor["97.5%", ])

colnames(pjor.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(pjor.summary, file = "Results/JAGS/pjor_summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.pjor <- matrix(NA, 195, 3)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.pjor) <- c("CI_2.5", "Median", "CI_97.5")
N.pjor
cbind(PJORmin, N.pjor)

##### DWRI #####
# Check fit
for(i in 1:3) bayesP.dwri4 <- mean(dwri_od[, "fit.new",][[i]] > dwri_od[, "fit",][[i]]) # 0.454 excellent
print(bayesP.dwri4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(dwri_od[, "fit",]), as.matrix(dwri_od[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=dwri_od[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), dig=3) # good convergence

Quants.dwri <- apply(as.matrix(dwri_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.dwri <- apply(as.matrix(dwri_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.dwri <- apply(as.matrix(dwri_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

Dwri.variables <- c("N-intercept", "Elevation", "TPI", "TWI", "Canopy", "Litter Depth", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

Dwri.summary <- data.frame(Dwri.variables, Means.dwri, SDs.dwri, Quants.dwri["2.5%", ], Quants.dwri["50%", ], Quants.dwri["97.5%", ])

colnames(Dwri.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(Dwri.summary, file = "Results/JAGS/Dwri_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.dwri <- matrix(NA, 195, 3)
for(i in 1:195){
  N.dwri[i, ] <- apply(as.matrix(dwri_od[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.dwri) <- c("CI_2.5", "Median", "CI_97.5")
N.dwri
cbind(DWRImin, N.dwri)


##### EWIL #####
# Check fit
for(i in 1:3) bayesP.ewil4 <- mean(ewil_od[, "fit.new",][[i]] > ewil_od[, "fit",][[i]]) # 0.454 excellent
print(bayesP.ewil4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(ewil_od[, "fit",]), as.matrix(ewil_od[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=ewil_od[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), dig=3) # good convergence

Quants.ewil <- apply(as.matrix(ewil_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.ewil <- apply(as.matrix(ewil_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.ewil <- apply(as.matrix(ewil_od[ , c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

ewil.variables <- c("N-intercept", "Elevation", "TPI", "TWI", "Canopy", "Litter Depth", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

ewil.summary <- data.frame(ewil.variables, Means.ewil, SDs.ewil, Quants.ewil["2.5%", ], Quants.ewil["50%", ], Quants.ewil["97.5%", ])

colnames(ewil.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(ewil.summary, file = "Results/JAGS/ewil_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.ewil <- matrix(NA, 195, 3)
for(i in 1:195){
  N.ewil[i, ] <- apply(as.matrix(ewil_od[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.ewil) <- c("CI_2.5", "Median", "CI_97.5")
N.ewil
cbind(EWILmin, N.ewil)
