
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

n.sites <- length(site.inits)

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = n.sites, 
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
test_nim <- nimbleCode({
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    # beta2.lam ~ dnorm(0, 0.01)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
      eps.lam[i] ~ dnorm(0, tau.site)
    }
    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    lambda[i] <- exp(alpha.lam + beta1.lam*elev[i] + eps.lam[i])
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    logit(p[i,j]) <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j]
    }
    }
    })

params <- c( "alpha.lam", 
             "beta1.lam", 
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta10.p",
             "sigma.site",
             "N",
             "p")

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       sigma.site = runif(1, 0.1, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta10.p = rnorm(1, 0.5, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

mcmc.output <- nimbleMCMC(test_nim, data = pjor.od.data, inits = inits,
                          monitors = params, thin = nt,
                          niter = ni+nb, nburnin = nb, nchains = 2,
                          summary = TRUE, WAIC = FALSE,
                          samplesAsCodaMCMC = TRUE)

library(coda)
pjor_od<- mcmc.list(mcmc.output$samples)
plot(pjor_od[ , c("alpha.lam", "beta1.lam", "alpha.p", "beta1.p", "beta2.p", "beta10.p", "sigma.site")]) # 
par(mfrow = c(1,1))

plot(mcmc.output$samples)


small_nim <- readBUGSmodel("Code/Jags_Models/small.txt", data = pjor.od.data, inits = NULL, dir = NULL, useInits = TRUE, debug = FALSE, returnComponents = FALSE, check = getNimbleOption("checkModel"), calculate = FALSE)

small_nim$initializeInfo()
# help("modelInitialization")

compileNimble(small_nim) 

samples <- nimbleMCMC(
  model = small_nim,
  monitors = params,
  niter = ni,
  nburnin = nb,
  thin = nt,
  nchains = 1)
