############################################
# Script for checking MCMC/HMC diagnostics
############################################

#----- Load Libraries ------

library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
# library(tidybayes)
library(loo)
# library(rstanarm)
library(bayesplot)
library(gridExtra)
library(stargazer)
library(readr)

#----- Set Defaults and Conditions ------

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

#----- Load Data and Results -----
load(file = "Data/Derived/settings.RData")
load("Data/Derived/stan_prep.RData")
site_od_full_pjor <- readRDS(file = "Results/Stan/final_od_pjor_hmc.Rds")

# load(file = "Results/JAGS/pjor2_mcmc_out.RData")

#----- Summarize coefficients -----

summary_table <- summary(site_od_full_pjor, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))$summary %>%
  data.frame() %>%
  # dplyr::mutate(parameter = rownames(.)) %>%
  select(-X25., -X75.)

write.csv(summary_table, file = "Results/Stan/pjor_summary.csv", row.names = TRUE)

# stargazer(summary(site_od_full_pjor)$summary,
#           type = "html",
#           out="Results/pjor_summary_table.doc")
          
          
#----- Plot parameter summaries -----

# ?`plot,stanfit-method`

# plot(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6"))

detection_labels <- c("Intercept", "Temperature", expression(Temperature^2), "Precip-24h", "Herbaceous", expression(Herbaceous^2), "Rel-Humidity")
N_labels <- c("Intercept", "Elevation", expression(Elevation^2), "TWI", "Litter", "Herbaceous", "Stream-Dist")

sims_mat <- as.matrix(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N")) 
 
mcmc_intervals(sims_mat, regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + rstan:::rstanvis_multiparam_theme()
mcmc_intervals(sims_mat, regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal()
mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal()

# stan_plot(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6")) + scale_y_discrete(labels = N_labels)




Quants.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = mean)

SDs.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = sd)

pjor.variables <- c("N-intercept", "Elevation", "Elevation^2", "TWI", "Litter Depth", "Ground Cover", "Stream Distance", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")

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


samples <- rstan::extract(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))
hist(samples$alpha0)

# Effect of Beta1 - elevation
elevation <- seq(450, 2025, length.out = 2000)
elevation_s <- (elevation - mean(Data$elev)) / sd(Data$elev)

N_hat_pjor <- matrix(NA, length(samples$alpha0), length(elevation_s))
for(i in 1:length(samples$alpha0)) {
  N_hat_pjor[i, ] <- exp(samples$alpha0[i] + samples$alpha1[i] * elevation_s + samples$alpha2[i] * elevation_s * elevation_s)
}
quants <- t(apply(as.matrix(N_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(elevation, quants)
colnames(quants_df) <- c("Elevation", "LCRI", "Median", "UCRI")

PJORmin <- apply(PJOR5, 1, max) 

N_pjor <- data.frame(summary(site_od_full_pjor, digits = 3, par = c("N"))$summary)
colnames(N_pjor) <- c("mean", "se_mean", "sd", "q2.5", "q25", "Median", "q75", "q97.5", "n_eff", "Rhat")

N_elev <- data.frame(Elevation = Data5$elev, N_pjor, PJORmin) # not sure if these line up the elevations correctly

gg_elev <- ggplot(quants_df, aes(Elevation, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance")  # + ggtitle(expression(paste(italic("Pethodon jordani")))) # 

gg_elev + geom_point(data = N_elev, aes(Elevation, Median))

ggplot(N_elev, aes(Elevation, Median)) + geom_crossbar(aes(ymin = q25, ymax = q75), width = 10, fill = "lightblue", colour = "lightblue") + geom_point() + ylab("Abundance") + geom_point(aes(Elevation, PJORmin, colour = "red"), size = 1) # + geom_legend("Max Count")











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

gg_gcover <- ggplot(quants_df, aes(gcover, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") + xlab("Ground cover (% area)") #+ coord_cartesian(ylim = c(0, 30)) # + geom_point(data = sna, aes(Elevation, Median))

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

## Ground cover detection curve

lp_hat_pjor <- matrix(NA, length(fit_pjor$alpha.p), length(gcover_s))
for(i in 1:length(fit_pjor$alpha.lam)) {
  lp_hat_pjor[i, ] <- fit_pjor$alpha.p[i] + fit_pjor$beta4.p[i] * gcover_s + fit_pjor$beta5.lam[i] * gcover_s * gcover_s
}
p_hat_pjor <- exp(lp_hat_pjor) / (1 + exp(lp_hat_pjor))
quants <- t(apply(as.matrix(p_hat_pjor), 2, FUN = stats::quantile, probs = c(0.025, 0.5, 0.975)))
quants_df <- data.frame(gcover, quants)
colnames(quants_df) <- c("gcover", "LCRI", "Median", "UCRI")

gg_gcover <- ggplot(quants_df, aes(gcover, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Detection probability") + xlab("Ground cover (% area)") #+ coord_cartesian(ylim = c(0, 30)) # + geom_point(data = sna, aes(Elevation, Median))


