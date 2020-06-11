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
library(ggpubr)

#----- Set Defaults and Conditions ------

if(!dir.exists("Results/Stan/Figures/")) dir.create("Results/Stan/Figures/", recursive = TRUE)

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

plot_cond <- function(fit, var = "elev", pars = c("alpha0", "alpha1", "alpha2"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95), link = "log") {
  samples <- rstan::extract(fit, pars = pars)
  lexpect <- matrix(NA_real_, length(samples[[1]]), length.out)
  if(is.null(range)) range <- base::range(as.matrix(Data[ , var]), na.rm = TRUE)
  x <- seq(range[1], range[2], length.out = length.out)
  xs <- (x - mean(as.matrix(Data[ , var]), na.rm = TRUE)) / sd(as.matrix(Data[ , var]), na.rm = TRUE)
  
  for(i in 1:length.out) {
    lexpect[ , i] <- samples[[pars[1]]] + samples[[pars[2]]] * xs[i] 
    if(length(pars) == 3) lexpect[ , i] <- lexpect[ , i] + samples[[pars[3]]] *  xs[i]^2
  }
  expect <- switch(link,
                   log = exp(lexpect),
                   logit = exp(lexpect) / (1 + exp(lexpect)))
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("X", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(X, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
}

plot_lin <- function(fit, var = "elev", pars = c("alpha0", "alpha6"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples[[1]]), length.out)
  if(is.null(range)) range <- base::range(Data[ , var], na.rm = TRUE)
  x <- seq(range[1], range[2], length.out = length.out)
  x_s <- (x - mean(as.matrix(Data[ , var]), na.rm = TRUE)) / sd(as.matrix(Data[ , var]), na.rm = TRUE)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * x_s[i])
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("X", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(X, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
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
 
# mcmc_intervals(sims_mat, regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + rstan:::rstanvis_multiparam_theme()
# mcmc_intervals(sims_mat, regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("P. jordani"))))
mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("P. jordani"))))
ggsave("Results/Stan/Figures/posteriors_p_pjor.pdf", dpi = 1000)

# stan_plot(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6")) + scale_y_discrete(labels = N_labels)

# mcmc_intervals(sims_mat, regex_pars = "alpha") + scale_y_discrete(labels = N_labels) + theme_bw_journal()
mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "alpha") + scale_y_discrete(labels = N_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("P. jordani"))))
ggsave("Results/Stan/Figures/posteriors_N_pjor.pdf", width = 8, height = 6, dpi = 1000)



# Quants.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))
# 
# Means.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = mean)
# 
# SDs.pjor <- apply(sims_mat[ , c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps")], 2, FUN = sd)
# 
# pjor.variables <- c("N-intercept", "Elevation", "Elevation^2", "TWI", "Litter Depth", "Ground Cover", "Stream Distance", "Site SD", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Detection SD", "Fit-Data", "Fit-Ideal")
# 
# pjor.summary <- data.frame(pjor.variables, Means.pjor, SDs.pjor, Quants.pjor["2.5%", ], Quants.pjor["50%", ], Quants.pjor["97.5%", ])
# 
# colnames(pjor.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")
# 
# write.table(pjor.summary, file = "Results/JAGS/pjor_summary.csv", sep = ",", col.names = NA, row.names = TRUE)
# 
# N.pjor <- matrix(NA, 159, 5)
# N.eff <- rep(NA, times = 159)
# for(i in 1:159){
#   # N.pjor[i, ] <- apply(sims_mat[, c(paste("N[", i, "]", sep = ""))], 1, FUN= quantile, probs = c(0.025, 0.1,  0.5, 0.9, 0.975))
#   N.pjor[i, ] <- quantile(sims_mat[, c(paste("N[", i, "]", sep = ""))], probs = c(0.025, 0.1,  0.5, 0.9, 0.975))
#   
#   # N.eff[i] <- effectiveSize(x=pjor_od2[ , c(paste("N[", i, "]", sep = ""))])
# }
# 
# colnames(N.pjor) <- c("CRI_2.5", "CRI_10", "Median", "CRI_90", "CRI_97.5")
# N.pjor
# cbind(PJORmin, N.pjor)







# combine chains into one for summarizing and plotting
# fit_pjor <- as.data.frame(do.call(rbind, pjor_od2))


samples <- rstan::extract(site_od_full_pjor, pars = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))
# hist(samples$alpha0)

# Effect of Beta1 - elevation

# s = extracted samples from rstan::extract
# pars = the parameter name or vector of parameter names  - need to include intercept - or jut go to design matrix setup - check pred for glm or lme4
# data = original data used for standardization

# fit = site_od_full_pjor
# data = Data$elev
# par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6")
# link = "log"
# data_new

# plot_eff <- function(fit, pars, data, data_new, link) {
#   if(!(link %in% c("log", "logit"))) stop("link only currently defined for log and logit")
#   
#   x <- range(data, na.rm = TRUE)
#   xs <- as.numeric(scale(data_new))
#   
#   samples <- rstan::extract(fit, par = par)
#   
#   lin_pred
#   
#   switch(link,
#          log = exp(lin_pred),
#          logit = exp(lin_pred) / (1 + exp(lin_pred)))
#   
# }

plot_elev <- function(fit, pars = c("alpha0", "alpha1", "alpha2"), data = Data, limits = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples$alpha0), length.out)
  if(is.null(limits)) limits <- base::range(Data$elev)
  x <- seq(limits[1], limits[2], length.out = length.out)
  xs <- (x - mean(Data$elev)) / sd(Data$elev)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * xs[i] + samples[[pars[3]]] *  xs[i]^2)
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("Elevation", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(Elevation, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + xlab("Elevation (m)") + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
return(gg)
}

plot_stream <- function(fit, pars = c("alpha0", "alpha6"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples$alpha0), length.out)
  if(is.null(range)) range <- base::range(Data$strm_dist)
  x <- seq(range[1], range[2], length.out = length.out)
  x_s <- (x - mean(Data$strm_dist)) / sd(Data$strm_dist)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * x_s[i])
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("x", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(x, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + xlab("Stream distance (m)") + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
}

plot_herb <- function(fit, pars = c("alpha0", "alpha5"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples$alpha0), length.out)
  if(is.null(range)) range <- base::range(Data$gcover, na.rm = TRUE)
  x <- seq(range[1], range[2], length.out = length.out)
  x_s <- (x - mean(Data$gcover, na.rm = TRUE)) / sd(Data$gcover, na.rm = TRUE)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * x_s[i])
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("x", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(x, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + xlab("Herbaceous ground cover (%)") + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
}

plot_litter <- function(fit, pars = c("alpha0", "alpha4"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples$alpha0), length.out)
  if(is.null(range)) range <- base::range(Data$litdepth)
  x <- seq(range[1], range[2], length.out = length.out)
  x_s <- (x - mean(Data$litdepth)) / sd(Data$litdepth)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * x_s[i])
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("x", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(x, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + xlab("Litter depth (mm)") + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
}

plot_twi <- function(fit, pars = c("alpha0", "alpha3"), data = Data, range = NULL, length.out = 1000, probs = c(0.05, 0.50, 0.95)) {
  samples <- rstan::extract(fit, pars = pars)
  expect <- matrix(NA_real_, length(samples$alpha0), length.out)
  if(is.null(range)) range <- base::range(Data$twi_10)
  x <- seq(range[1], range[2], length.out = length.out)
  x_s <- (x - mean(Data$twi_10)) / sd(Data$twi_10)
  
  for(i in 1:length.out) {
    expect[ , i] <- exp(samples[[pars[1]]] + samples[[pars[2]]] * x_s[i])
  }
  quants <- t(apply(expect, 2, FUN = stats::quantile, probs = probs))
  quants_df <- data.frame(x, quants)
  colnames(quants_df) <- c("x", "LCRI", "Median", "UCRI")
  
  gg <- ggplot(quants_df, aes(x, Median)) + geom_line() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha=0.3) + xlab("TWI") + ylab("Abundance") # + ggtitle(expression(paste(italic("Pethodon jordani")))) #
  return(gg)
}

gg_N_elev <- plot_elev(site_od_full_pjor) + coord_cartesian(ylim = c(0, 115))
gg_N_twi <- plot_twi(site_od_full_pjor) + coord_cartesian(ylim = c(0, 115))
gg_N_litter <- plot_litter(site_od_full_pjor) + coord_cartesian(ylim = c(0, 115))
gg_N_herb <- plot_herb(site_od_full_pjor) + coord_cartesian(ylim = c(0, 115))
gg_N_stream <- plot_stream(site_od_full_pjor) + coord_cartesian(ylim = c(0, 115))

library(ggpubr)
gg_N <- ggarrange(gg_N_elev, gg_N_litter + ylab(""), gg_N_herb, gg_N_stream + ylab(""), 
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
gg_N <- annotate_figure(gg_N, top = text_grob(expression(paste(italic("P. jordani"))), face = "bold", size = 15))
gg_N
ggsave(filename = "Results/Stan/Figures/abundance_pjor.pdf", width = 8, height = 6, gg_N, dpi = 1000)

# same thing for detection
gg_p_temp <- plot_cond(site_od_full_pjor, var = c("temp1", "temp2", "temp3", "temp4", "temp5", "temp6"), pars = c("beta0", "beta1", "beta2"), link = "logit") + xlab("Temperature (C)") + ylab("Prob. of detection")
gg_p_precip <- plot_cond(site_od_full_pjor, var = c("precip1", "precip2", "precip3", "precip4", "precip5", "precip6"), pars = c("beta0", "beta3"), link = "logit") + xlab("Precip. in past 24-hr (mm)") + ylab("Prob. of detection")
gg_p_ground <- plot_cond(site_od_full_pjor, var = c("ground1", "ground2", "ground3", "ground4", "ground5", "ground6"), pars = c("beta0", "beta4", "beta5"), link = "logit") + xlab("Herbaceous ground cover (%)") + ylab("Prob. of detection")
gg_p_RH <- plot_cond(site_od_full_pjor, var = c("RH1", "RH2", "RH3", "RH4", "RH5", "RH6"), pars = c("beta0", "beta6"), link = "logit") + xlab("Relative humidity (%)") + ylab("Prob. of detection")

gg_p <- ggarrange(gg_p_temp + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_precip + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab")
                  gg_p_ground + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_RH + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab") 
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
gg_p <- annotate_figure(gg_p, top = text_grob(expression(paste(italic("P. jordani"))), face = "bold", size = 15)) # rel(1.5))) # adjust sizes
gg_p
ggsave(plot = gg_p, filename = "Results/Stan/Figures/detection_pjor.png", width = 8, height = 6)
ggsave(plot = gg_p, filename = "Results/Stan/Figures/detection_pjor.pdf", width = 8, height = 6, dpi = 1000)

#----- EWIL -----

site_od_full_ewil <- readRDS(file = "Results/Stan/final_od_ewil_hmc.Rds")

#----- Summarize coefficients -----

summary_table <- summary(site_od_full_ewil, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))$summary %>%
  data.frame() %>%
  # dplyr::mutate(parameter = rownames(.)) %>%
  select(-X25., -X75.)

write.csv(summary_table, file = "Results/Stan/ewil_summary.csv", row.names = TRUE)

#----- Plot parameter summaries -----

sims_mat <- as.matrix(site_od_full_ewil, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N")) 

mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("E. wilderae"))))
ggsave("Results/Stan/Figures/posteriors_p_ewil.pdf", dpi = 1000)

mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "alpha") + scale_y_discrete(labels = N_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("E. wilderae"))))
ggsave("Results/Stan/Figures/posteriors_N_ewil.pdf", dpi = 1000)

samples <- rstan::extract(site_od_full_ewil, pars = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))

gg_N_elev <- plot_elev(site_od_full_ewil) + coord_cartesian(ylim = c(0, 20))
gg_N_twi <- plot_twi(site_od_full_ewil) + coord_cartesian(ylim = c(0, 20))
gg_N_litter <- plot_litter(site_od_full_ewil) + coord_cartesian(ylim = c(0, 20))
gg_N_herb <- plot_herb(site_od_full_ewil) + coord_cartesian(ylim = c(0, 20))
gg_N_stream <- plot_stream(site_od_full_ewil) + coord_cartesian(ylim = c(0, 20))

gg_N <- ggarrange(gg_N_elev, gg_N_litter + ylab(""), gg_N_herb, gg_N_stream + ylab(""), 
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
gg_N <- annotate_figure(gg_N, top = text_grob(expression(paste(italic("E. wilderae"))), face = "bold", size = 15))
gg_N
ggsave(plot = gg_N, filename = "Results/Stan/Figures/abundance_ewil.pdf", width = 8, height = 6, dpi = 1000)
rm(gg_N)

# same thing for detection
gg_p_temp <- plot_cond(site_od_full_ewil, var = c("temp1", "temp2", "temp3", "temp4", "temp5", "temp6"), pars = c("beta0", "beta1", "beta2"), link = "logit") + xlab("Temperature (C)") + ylab("Prob. of detection")
gg_p_precip <- plot_cond(site_od_full_ewil, var = c("precip1", "precip2", "precip3", "precip4", "precip5", "precip6"), pars = c("beta0", "beta3"), link = "logit") + xlab("Precip. in past 24-hr (mm)") + ylab("Prob. of detection")
gg_p_ground <- plot_cond(site_od_full_ewil, var = c("ground1", "ground2", "ground3", "ground4", "ground5", "ground6"), pars = c("beta0", "beta4", "beta5"), link = "logit") + xlab("Herbaceous ground cover (%)") + ylab("Prob. of detection")
gg_p_RH <- plot_cond(site_od_full_ewil, var = c("RH1", "RH2", "RH3", "RH4", "RH5", "RH6"), pars = c("beta0", "beta6"), link = "logit") + xlab("Relative humidity (%)") + ylab("Prob. of detection")

gg_p <- ggarrange(gg_p_temp + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_precip + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab")
                  gg_p_ground + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_RH + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab") 
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
gg_p <- annotate_figure(gg_p, top = text_grob(expression(paste(italic("E. wilderae"))), face = "bold", size = 15))
gg_p
ggsave(plot = gg_p, filename = "Results/Stan/Figures/detection_ewil.pdf", width = 8, height = 6, dpi = 1000)
rm(gg_p)

#----- DWRI -----

site_od_full_dwri <- readRDS(file = "Results/Stan/final_od_dwri_hmc.Rds")

#----- Summarize coefficients -----

summary_table <- summary(site_od_full_dwri, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))$summary %>%
  data.frame() %>%
  # dplyr::mutate(parameter = rownames(.)) %>%
  select(-X25., -X75.)

write.csv(summary_table, file = "Results/Stan/dwri_summary.csv", row.names = TRUE)

#----- Plot parameter summaries -----

sims_mat <- as.matrix(site_od_full_dwri, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N")) 

mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "beta") + scale_y_discrete(labels = detection_labels) + theme_bw_journal() + ggtitle(expression(paste(italic("D. wrighti"))))
ggsave("Results/Stan/Figures/posteriors_p_dwri.pdf", dpi = 1000)

mcmc_areas(sims_mat, area_method = "scaled height", prob = 0.5, prob_outer = 0.95, point_est = "median", regex_pars = "alpha") + scale_y_discrete(labels = N_labels) + ggtitle(expression(paste(italic("D. wrighti")))) + theme_bw_journal()
ggsave("Results/Stan/Figures/posteriors_N_dwri.pdf", dpi = 1000)

samples <- rstan::extract(site_od_full_dwri, pars = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p", "sd_eps", "N"))

gg_N_elev <- plot_elev(site_od_full_dwri) + coord_cartesian(ylim = c(0, 20))
gg_N_twi <- plot_twi(site_od_full_dwri) + coord_cartesian(ylim = c(0, 20))
gg_N_litter <- plot_litter(site_od_full_dwri) + coord_cartesian(ylim = c(0, 20))
gg_N_herb <- plot_herb(site_od_full_dwri) + coord_cartesian(ylim = c(0, 20))
gg_N_stream <- plot_stream(site_od_full_dwri) + coord_cartesian(ylim = c(0, 20))

gg_N <- ggarrange(gg_N_elev, gg_N_litter + ylab(""), gg_N_herb, gg_N_stream + ylab(""), 
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
gg_N <- annotate_figure(gg_N, top = text_grob(expression(paste(italic("D. wrighti"))), face = "bold", size = 15))
gg_N
ggsave(plot = gg_N, filename = "Results/Stan/Figures/abundance_dwri.pdf", width = 8, height = 6, dpi = 1000)
rm(gg_N)

# same thing for detection
gg_p_temp <- plot_cond(site_od_full_dwri, var = c("temp1", "temp2", "temp3", "temp4", "temp5", "temp6"), pars = c("beta0", "beta1", "beta2"), link = "logit") + xlab("Temperature (C)") + ylab("Prob. of detection")
gg_p_precip <- plot_cond(site_od_full_dwri, var = c("precip1", "precip2", "precip3", "precip4", "precip5", "precip6"), pars = c("beta0", "beta3"), link = "logit") + xlab("Precip. in past 24-hr (mm)") + ylab("Prob. of detection")
gg_p_ground <- plot_cond(site_od_full_dwri, var = c("ground1", "ground2", "ground3", "ground4", "ground5", "ground6"), pars = c("beta0", "beta4", "beta5"), link = "logit") + xlab("Herbaceous ground cover (%)") + ylab("Prob. of detection")
gg_p_RH <- plot_cond(site_od_full_dwri, var = c("RH1", "RH2", "RH3", "RH4", "RH5", "RH6"), pars = c("beta0", "beta6"), link = "logit") + xlab("Relative humidity (%)") + ylab("Prob. of detection")

gg_p <- ggarrange(gg_p_temp + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_precip + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab")
                  gg_p_ground + coord_cartesian(ylim = c(0, 0.6)), 
                  gg_p_RH + coord_cartesian(ylim = c(0, 0.6)) + ylab(""), # + rremove("y.text") + rremove("ylab") 
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
gg_p <- annotate_figure(gg_p, top = text_grob(expression(paste(italic("D. wrighti"))), face = "bold", size = 15))
gg_p
ggsave(plot = gg_p, filename = "Results/Stan/Figures/detection_dwri.pdf", width = 8, height = 6, dpi = 1000)


#---------- Combined table ------------

detection_labels <- c("Intercept", "Temperature", "Temperature^2", "Precip-24h", "Herbaceous", "Herbaceous^2", "Rel-Humidity", "Random Obs. SD")

N_labels <- c("Intercept", "Elevation", "Elevation^2", "TWI", "Litter", "Herbaceous", "Stream-Dist", "Site SD")

library(tibble)
pjor_sum <- summary(site_od_full_pjor, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "sd_eps", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p"))$summary %>%
  data.frame() %>%
  dplyr::mutate(parameter = rownames(.)) %>%
  add_column(Variable = c(N_labels, detection_labels)) %>%
  select(Variable, Mean = mean, X2.5., X97.5.)

dwri_sum <- summary(site_od_full_dwri, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "sd_eps", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p"))$summary %>%
  data.frame() %>%
  dplyr::mutate(parameter = rownames(.)) %>%
  add_column(Variable = c(N_labels, detection_labels)) %>%
  select(Variable, Mean = mean, X2.5., X97.5.)

ewil_sum <- summary(site_od_full_ewil, digits = 3, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "sd_eps", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p"))$summary %>%
  data.frame() %>%
  dplyr::mutate(parameter = rownames(.)) %>%
  add_column(Variable = c(N_labels, detection_labels)) %>%
  select(Variable, Mean = mean, X2.5., X97.5.)

summary_table_3sp <- data.frame(pjor_sum, dwri_sum[2:4], ewil_sum[2:4])
write.csv(summary_table_3sp, file = "Results/Stan/summary_table_3sp.csv", row.names = FALSE)

# sims_mat <- as.matrix(site_od_full_pjor, par = c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "sd_eps", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sd_p"))


