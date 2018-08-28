# Summarize JAGS results

# Load libraries
library(rjags)
library(coda)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# Set defaults
theme_bw_journal <- function (base_family = "") {
  theme_grey(base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = rel(1.2)), 
      axis.title = element_text(size =rel(1.25)), 
      axis.ticks = element_line(colour = "black"),
      legend.key = element_rect(colour = "grey80"),
      panel.background = element_rect(fill = "white",
                                      colour = NA), 
      panel.border = element_rect(fill = NA, 
                                  colour = "grey50"), 
      panel.grid.major = element_line(colour = "grey90", 
                                      size = 0.2), 
      panel.grid.minor = element_line(colour = "grey98",
                                      size = 0.5), 
      strip.background = element_rect(fill = "grey80", 
                                      colour = "grey50", 
                                      size = 0.),
      plot.title = element_text(hjust = 0.5,
                                size = rel(1.5)))
  }
theme_set(theme_bw_journal())
# theme_update(plot.title = element_text(hjust = 0.5))

# Load Data
load("Data/Processed/jags_prep.RData")
load(file = "Results/JAGS/pjor_mcmc_out.RData")
load(file = "Results/JAGS/pjor2_mcmc_out.RData")
load(file = "Results/JAGS/dwri_mcmc_out.RData")
load(file = "Results/JAGS/dwri2_mcmc_out.RData")
load(file = "Results/JAGS/ewil_mcmc_out.RData")


##### PJOR #####
# Check fit
for(i in 1:3) bayesP.pjor4 <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # ~0.4 good
print(bayesP.pjor4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # 
abline(0, 1, col = 'red')

plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")])

print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")]), dig=3) # good convergence

print(effectiveSize(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[10]", "N[50]", "N[75]", "N[100]", "N[125]", "N[150]", "N[175]", "N[195]")]), dig=3)

Quants.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.pjor <- apply(as.matrix(pjor_od[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

pjor.variables <- c("N-intercept", "Elevation", "TWI", "Litter Depth", "Ground Cover", "Stream Distance", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

pjor.summary <- data.frame(pjor.variables, Means.pjor, SDs.pjor, Quants.pjor["2.5%", ], Quants.pjor["50%", ], Quants.pjor["97.5%", ])

colnames(pjor.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(pjor.summary, file = "Results/JAGS/pjor_summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.pjor <- matrix(NA, 195, 3)
N.eff <- rep(NA, times = 195)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))

N.eff[i] <- effectiveSize(x=pjor_od[ , c(paste("N[", i, "]", sep = ""))])
}

colnames(N.pjor) <- c("CI_2.5", "Median", "CI_97.5")
N.pjor
cbind(PJORmin, N.pjor)

# combine chains into one for summarizing and plotting
fit_pjor <- as.data.frame(do.call(rbind, pjor_od))

hist(fit_pjor$alpha.lam)

# Effect of Beta1 - elevation
elevation <- seq(450, 2025, length.out = 2000)
elevation_s <- (elevation - mean(Data$elev)) / sd(Data$elev)

N_hat_pjor <- matrix(NA, length(fit_pjor$beta1.lam), length(elevation_s))
for(i in 1:length(fit_pjor$beta1.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta1.lam[i] * elevation_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(elevation, quants)
colnames(quants_df) <- c("Elevation", "LCRI", "Median", "UCRI")

# sna <- data.frame(Elevation = Data$elev, N.pjor) # not sure if these line up the elevations correctly

gg_elev <- ggplot(quants_df, aes(Elevation, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) # + geom_point(data = sna, aes(Elevation, Median))

# Effect of Beta4 - litter depth
depth <- seq(min(Data$litdepth), max(Data$litdepth), length.out = 2000)
depth_s <- (depth - mean(Data$litdepth)) / sd(Data$litdepth)

N_hat_pjor <- matrix(NA, length(fit_pjor$alpha.lam), length(depth_s))
for(i in 1:length(fit_pjor$alpha.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta4.lam[i] * depth_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(depth, quants)
colnames(quants_df) <- c("Depth", "LCRI", "Median", "UCRI")

gg_depth <- ggplot(quants_df, aes(Depth, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") + xlab("Litter depth (mm)") + coord_cartesian(ylim = c(0, 30)) # + geom_point(data = sna, aes(Elevation, Median))

# Effect of Beta5 - ground cover
gcover <- seq(min(Data$gcover, na.rm = TRUE), max(Data$gcover, na.rm = TRUE), length.out = 2000)
gcover_s <- (gcover - mean(Data$gcover, na.rm = TRUE)) / sd(Data$gcover, na.rm = TRUE)

N_hat_pjor <- matrix(NA, length(fit_pjor$alpha.lam), length(gcover_s))
for(i in 1:length(fit_pjor$alpha.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta5.lam[i] * gcover_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(gcover, quants)
colnames(quants_df) <- c("gcover", "LCRI", "Median", "UCRI")

gg_gcover <- ggplot(quants_df, aes(gcover, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") + xlab("Ground Gover (% area)") #+ coord_cartesian(ylim = c(0, 30)) # + geom_point(data = sna, aes(Elevation, Median))

pjor_N_grid <- arrangeGrob(gg_elev + coord_cartesian(ylim = c(0, 1000)), 
                      gg_depth + coord_cartesian(ylim = c(0, 1000)) + ylab(""), 
                      gg_gcover + coord_cartesian(ylim = c(0, 1000)) + ylab(""), 
                      ncol = 3, 
                      nrow = 1, 
                      widths = c(3.5, 3.5, 3.5)) #, 
                      # heights = c(rep(5, times=5)))

plot(pjor_N_grid)

if(!dir.exists("Results/JAGS/Figures")) dir.create("Results/JAGS/Figures", recursive = TRUE)

ggsave(file = "Results/JAGS/Figures/pjor_N_grid.pdf", pjor_N_grid, dpi = 1000)

##### PJOR Elev2 #####
# Check fit
for(i in 1:3) bayesP.pjor4 <- mean(pjor_od2[, "fit.new",][[i]] > pjor_od2[, "fit",][[i]]) # ~0.4 good
print(bayesP.pjor4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od2[, "fit",]), as.matrix(pjor_od2[, "fit.new",])) # 
abline(0, 1, col = 'red')

plot(pjor_od2[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")]) # 
par(mfrow = c(1,1))
summary(pjor_od2[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")])

print(gelman.diag(x=pjor_od2[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")]), dig=3) # good convergence

print(effectiveSize(x=pjor_od2[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new", "N[1]", "N[10]", "N[50]", "N[75]", "N[100]", "N[125]", "N[150]", "N[175]", "N[195]")]), dig=3)

Quants.pjor <- apply(as.matrix(pjor_od2[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.pjor <- apply(as.matrix(pjor_od2[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.pjor <- apply(as.matrix(pjor_od2[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "sigma.site", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

pjor.variables <- c("N-intercept", "Elevation", "TWI", "Litter Depth", "Ground Cover", "Stream Distance", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

pjor.summary <- data.frame(pjor.variables, Means.pjor, SDs.pjor, Quants.pjor["2.5%", ], Quants.pjor["50%", ], Quants.pjor["97.5%", ])

colnames(pjor.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(pjor.summary, file = "Results/JAGS/pjor_summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.pjor <- matrix(NA, 195, 5)
N.eff <- rep(NA, times = 195)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od2[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.1,  0.5, 0.9, 0.975))
  
  N.eff[i] <- effectiveSize(x=pjor_od2[ , c(paste("N[", i, "]", sep = ""))])
}

colnames(N.pjor) <- c("CRI_2.5", "CRI_10", "Median", "CRI_90", "CRI_97.5")
N.pjor
cbind(PJORmin, N.pjor)

# combine chains into one for summarizing and plotting
fit_pjor <- as.data.frame(do.call(rbind, pjor_od2))

hist(fit_pjor$alpha.lam)

# Effect of Beta1 - elevation
elevation <- seq(450, 2025, length.out = 2000)
elevation_s <- (elevation - mean(Data$elev)) / sd(Data$elev)

N_hat_pjor <- matrix(NA, length(fit_pjor$beta1.lam), length(elevation_s))
for(i in 1:length(fit_pjor$beta1.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta1.lam[i] * elevation_s + fit_pjor$beta2.lam[i] * elevation_s * elevation_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(elevation, quants)
colnames(quants_df) <- c("Elevation", "LCRI", "Median", "UCRI")

N_elev <- data.frame(Elevation = Data$elev, N.pjor, PJORmin) # not sure if these line up the elevations correctly

gg_elev <- ggplot(quants_df, aes(Elevation, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance")  # + ggtitle(expression(paste(italic("Pethodon jordani")))) # 

gg_elev + geom_point(data = N_elev, aes(Elevation, Median))

ggplot(N_elev, aes(Elevation, Median)) + geom_crossbar(aes(ymin = CRI_10, ymax = CRI_90), width = 10, fill = "lightblue", colour = "lightblue") + geom_point() + ylab("Abundance") + geom_point(aes(Elevation, PJORmin, colour = "red"), size = 1) # + geom_legend("Max Count")

# Effect of Beta4 - litter depth
depth <- seq(min(Data$litdepth), max(Data$litdepth), length.out = 2000)
depth_s <- (depth - mean(Data$litdepth)) / sd(Data$litdepth)

N_hat_pjor <- matrix(NA, length(fit_pjor$alpha.lam), length(depth_s))
for(i in 1:length(fit_pjor$alpha.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta4.lam[i] * depth_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(depth, quants)
colnames(quants_df) <- c("Depth", "LCRI", "Median", "UCRI")

gg_depth <- ggplot(quants_df, aes(Depth, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") + xlab("Litter depth (mm)") + coord_cartesian(ylim = c(0, 125)) # + geom_point(data = sna, aes(Elevation, Median))

# Effect of Beta5 - ground cover
gcover <- seq(min(Data$gcover, na.rm = TRUE), max(Data$gcover, na.rm = TRUE), length.out = 2000)
gcover_s <- (gcover - mean(Data$gcover, na.rm = TRUE)) / sd(Data$gcover, na.rm = TRUE)

N_hat_pjor <- matrix(NA, length(fit_pjor$alpha.lam), length(gcover_s))
for(i in 1:length(fit_pjor$alpha.lam)) {
  N_hat_pjor[i, ] <- exp(fit_pjor$alpha.lam[i] + fit_pjor$beta5.lam[i] * gcover_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(gcover, quants)
colnames(quants_df) <- c("gcover", "LCRI", "Median", "UCRI")

gg_gcover <- ggplot(quants_df, aes(gcover, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") + xlab("Ground Gover (% area)") #+ coord_cartesian(ylim = c(0, 30)) # + geom_point(data = sna, aes(Elevation, Median))

pjor2_N_grid <- arrangeGrob(gg_elev + coord_cartesian(ylim = c(0, 250)), 
                           gg_depth + coord_cartesian(ylim = c(0, 250)) + ylab(""), 
                           gg_gcover + coord_cartesian(ylim = c(0, 250)) + ylab(""), 
                           ncol = 3, 
                           nrow = 1, 
                           widths = c(3.5, 3.5, 3.5)) #, 
# heights = c(rep(5, times=5)))

plot(pjor2_N_grid)

if(!dir.exists("Results/JAGS/Figures")) dir.create("Results/JAGS/Figures", recursive = TRUE)

ggsave(file = "Results/JAGS/Figures/pjor_N_grid.pdf", pjor2_N_grid, dpi = 1000)

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

write.table(Dwri.summary, file = "Results/JAGS/Dwri_summary.csv", sep = ",", col.names = NA, row.names = TRUE)

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

write.table(ewil.summary, file = "Results/JAGS/ewil_summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.ewil <- matrix(NA, 195, 3)
for(i in 1:195){
  N.ewil[i, ] <- apply(as.matrix(ewil_od[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.ewil) <- c("CI_2.5", "Median", "CI_97.5")
N.ewil
cbind(EWILmin, N.ewil)
