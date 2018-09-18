
library(rjags)
library(coda)
library(parallel)

load("Data/Processed/jags_prep.RData")

tpi <- (Data$tpi_100 - mean(Data$tpi_100))/sd(Data$tpi_100)

stream <- (Data$strm_dist - mean(Data$strm_dist)) / sd(Data$strm_dist)

testing <- FALSE
if(testing) {
  nb = 5000
  ni = 1000
  nt = 1
} else {
  nb = 100000
  ni = 400000
  nt = 160
}

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
                     elev = elev, 
                     elev2 = elev*elev, 
                     slope = slope,
                     slope2 = slope*slope, 
                     aspectN = aspectN,
                     aspectE = aspectE,
                     twi = twi,
                     tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover*gcover,
                     litterdepth = litterdepth,
                     stream = stream,
                     streamelev = stream*elev,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
                     elev = elev, 
                     elev2 = elev*elev, 
                     slope = slope,
                     slope2 = slope*slope, 
                     aspectN = aspectN,
                     aspectE = aspectE,
                     twi = twi,
                     tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover*gcover,
                     litterdepth = litterdepth,
                     stream = stream,
                     streamelev = stream*elev,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

ewil.od.data <- list(C = as.matrix(EWIL[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
                     elev = elev, 
                     elev2 = elev*elev, 
                     slope = slope,
                     slope2 = slope*slope, 
                     streamelev = stream*elev,
                     aspectN = aspectN,
                     aspectE = aspectE,
                     twi = twi,
                     tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover*gcover,
                     litterdepth = litterdepth,
                     stream = stream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])


####### PJOR ######

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             #"beta2.lam", 
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12591)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od <- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")])

save(pjor_od, file = "Results/JAGS/pjor_mcmc_out.RData")
# saveRDS(pjor_od, file = "Results/JAGS/pjor_mcmc_out.rds")

####### PJOR Elev^2 ######

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12592)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od2 <- mcmc.list(out)
plot(pjor_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new", "N[1]", "N[60]", "N[100]", "N[150]")])

save(pjor_od2, file = "Results/JAGS/pjor2_mcmc_out.RData")
# saveRDS(pjor_od2, file = "Results/JAGS/pjor2_mcmc_out.rds")

####### PJOR Stream-Elevation ######

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             #"beta2.lam", 
             "beta6.lam",
             "beta8.lam",
             "beta9.lam",
             "beta11.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12591)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/stream-elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od_stream <- mcmc.list(out)
plot(pjor_od_stream[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od_stream[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")])

save(pjor_od_stream, file = "Results/JAGS/pjor_stream_mcmc_out.RData")
# saveRDS(pjor_od, file = "Results/JAGS/pjor_mcmc_out.rds")

####### PJOR Stream-Elevation ######

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta6.lam",
             "beta8.lam",
             "beta9.lam",
             "beta11.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12591)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/stream-elev2_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od_stream2 <- mcmc.list(out)
plot(pjor_od_stream2[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od_stream2[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta8.lam", "beta9.lam", "beta11.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new", "N[1]", "N[150]")])

save(pjor_od_stream, file = "Results/JAGS/pjor_stream2_mcmc_out.RData")
# saveRDS(pjor_od, file = "Results/JAGS/pjor_mcmc_out.rds")

##### DWRI #####

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             #"beta2.lam", 
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12593)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_elev_od.txt", dwri.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
dwri_od <- mcmc.list(out)
plot(dwri_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(dwri_od, file = "Results/JAGS/dwri_mcmc_out.RData")
# saveRDS(dwri_od, file = "Results/JAGS/dwri_mcmc_out.rds")

##### DWRI Elev^2 #####

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12590)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_od.txt", dwri.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
dwri_od2 <- mcmc.list(out)
plot(dwri_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(dwri_od2, file = "Results/JAGS/dwri2_mcmc_out.RData")
# saveRDS(dwri_od2, file = "Results/JAGS/dwri2_mcmc_out.rds")


##### EWIL #####

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 6, 1),
       beta2.lam = rnorm(1, -6, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12590)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_od.txt", ewil.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
ewil_od <- mcmc.list(out)
plot(ewil_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(ewil_od, file = "Results/JAGS/ewil_mcmc_out.RData")
# saveRDS(ewil_od, file = "Results/JAGS/ewil_mcmc_out.rds")

rm(out)
gc()


##### EWIL uniform priors #####

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = runif(1, 4, 7),
       beta2.lam = rnorm(1, -8, 4),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12590)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_od_unif.txt", ewil.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
ewil_od_unif <- mcmc.list(out)
plot(ewil_od_unif[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_od_unif[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(ewil_od_unif, file = "Results/JAGS/ewil_unif_mcmc_out.RData")
# saveRDS(ewil_od_unif, file = "Results/JAGS/ewil_mcmc_out.rds")

rm(out)
gc()

##### EWIL Stream^2 #####

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam",
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "beta7.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 12590)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_stream2_od.txt", ewil.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
ewil_od2 <- mcmc.list(out)
plot(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(ewil_od2, file = "Results/JAGS/ewil_stream2_mcmc_out.RData")
# saveRDS(ewil_od, file = "Results/JAGS/ewil_mcmc_out.rds")

##### EWIL Full #####

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             "beta10.lam",
             "beta11.lam",
             "beta12.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             "eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

ewil.od.data <- list(C = as.matrix(EWIL[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
                     elev = elev, 
                     elev2 = elev*elev, 
                     slope = slope,
                     slope2 = slope*slope, 
                     aspectN = aspectN,
                     aspectE = aspectE,
                     ltwi = twi,
                     tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover*gcover,
                     litterdepth = litterdepth,
                     lstream = stream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_od.txt", ewil.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
ewil_od_full <- mcmc.list(out)
plot(ewil_od_full[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_od_full[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

save(ewil_od_full, file = "Results/JAGS/ewil_full_mcmc_out.RData")
# saveRDS(ewil_od, file = "Results/JAGS/ewil_mcmc_out.rds")

rm(out)
gc()





save(pjor_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/pjor_mcmc_out.RData")
saveRDS(pjor_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/pjor_mcmc_out.rds")

save(dwri_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/dwri_mcmc_out.RData")
saveRDS(dwri_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/dwri_mcmc_out.rds")

save(ewil_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/ewil_mcmc_out.RData")
saveRDS(ewil_od, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/ewil_mcmc_out.rds")

save(ewil_od2, file = "/Users/djhocking/OneDrive\ -\ Frostburg\ State\ University/Elevation_2012_Results/ewil_stream2_mcmc_out.RData")
