
# install.packages("nimble", type = "source", repos = "http://r-nimble.org")
library(nimble)
library(coda)

nimbleOptions(showCompilerOutput = TRUE) 

load("Data/Processed/jags_prep.RData")

testing <- TRUE
if(testing) {
  nb = 5000
  ni = 1000
  nt = 1
} else {
  nb = 50000
  ni = 100000
  nt = 50
}

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     aspectN = aspectN,
                     aspectE = aspectE,
                     ltwi = ltwi,
                     tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover2,
                     litterdepth = litterdepth,
                     lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])



#----------min elev2 model site abund no od in detection----Doesn't Fit--------
# random effect of site on abundance in detection
# Define model
sink("Code/Jags_Models/small.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    lambda[i] <- exp(beta1.lam*elev[i] + beta2.lam*elev2[i] + eps.lam[site[i]])
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    logit(p[i,j]) <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j]
    }
    }
}
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0.1, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta10.p",
             "eps.lam",
             # "delta.p",
             "sigma.site",
             # "sigma.p",
             "N",
             "p")

small_nim <- readBUGSmodel("Code/Jags_Models/small.txt", data = pjor.od.data, inits = inits(), dir = NULL, useInits = TRUE, debug = FALSE, returnComponents = FALSE, check = getNimbleOption("checkModel"), calculate = FALSE)

small_nim$initializeInfo()
# help("modelInitialization")

# compileNimble(small_nim) 

samples <- nimbleMCMC(
  model = small_nim,
  monitors = params,
  niter = ni+nb,
  nburnin = nb,
  thin = nt,
  nchains = 3,
  WAIC = FALSE,
  samplesAsCodaMCMC = TRUE)

# Results
pjor_od<- mcmc.list(samples)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta10.p", "sigma.site")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta10.p", "sigma.site")])








#----------Full model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
sink("Code/Jags_Models/full_elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)#T(-3, 5)
   
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta10.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    beta12.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- 1/(sigma.p*sigma.p)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    lambda[i] <- exp(beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]])
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    logit(p[i,j]) <- beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]

    }
    }
    
    }
    ", fill = TRUE)
sink()


params <- c( "alpha.lam", 
             "beta1.lam", 
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
             "p")

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
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


full_od_nim <- readBUGSmodel("Code/Jags_Models/full_elev_od.txt", data = pjor.od.data, inits = inits(), dir = NULL, useInits = TRUE, debug = FALSE, returnComponents = FALSE, check = getNimbleOption("checkModel"), calculate = FALSE)

full_od_nim$initializeInfo()

samples <- nimbleMCMC(
  model = full_od_nim,
  monitors = params,
  niter = ni + nb,
  nburnin = nb,
  thin = nt,
  nchains = 3,
  WAIC = FALSE,
  samplesAsCodaMCMC = TRUE)

# Results
pjor_od <- mcmc.list(samples)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "alpha.p", "beta1.p", "beta2.p", "beta10.p", "sigma.site")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta10.p", "sigma.site")])






#----------Overdispersion in detection simple spatial model no elevation^2-------------

params <- c( "alpha.lam", 
             "beta1.lam", 
             # "beta2.lam", 
             "beta3.lam",
             # "beta4.lam",
             # "beta5.lam",
             # "beta6.lam",
             # "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             # "beta10.lam",
             # "beta11.lam",
             # "beta12.lam",
             # "beta13.lam",
             "theta",
             "w",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             # "eps.lam",
             "delta.p",
             # "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")


elev_sp_od_nim <- readBUGSmodel("Code/Jags_Models/elev_sp_od.txt", data = pjor.od.data, inits = NULL, dir = NULL,
              useInits = TRUE, debug = FALSE, returnComponents = FALSE,
              check = getNimbleOption("checkModel"), calculate = FALSE)

compileNimble(elev_sp_od_nim) 

samples <- nimbleMCMC(
  model = elev_sp_od_nim,
  monitors = params,
  niter = ni,
  nburnin = nb,
  thin = nt)

