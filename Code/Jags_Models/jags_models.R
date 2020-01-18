
if(!dir.exists("Code/Jags_Models")) dir.create("Code/Jags_Models")

#----------Overdispersion in both abund and detection + kriging + no covariates------------
# Define model
sink("Code/Jags_Models/sp_od_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
  
    
    alpha.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[], Omega[,]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    # theta ~ dgamma(0.01, 0.01) # spatial correlation rate
# theta ~ dunif(0.0001, 1000) # spatial correlation rate
    theta ~ dnorm(0, 0.01)T(0.00001, 1000) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
    mu.lam[i] ~ dnorm(alpha.lam, tau.lam)
    # muW[i] <- mu.lam[i]
muW[i] <- 0
    }
    
    tau.lam <- 1 / (sigma.lam * sigma.lam)
    sigma.lam ~ dunif(0, 3)
    
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- mu.lam[i] + W[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       # beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       # beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       theta = runif(1, 2, 50),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "theta",
             "w",
             "alpha.p", 
             # "eps.lam",
             "sigma.lam",
             "delta.p",
             # "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     # elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover * gcover,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/sp_od_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "alpha.p", "theta", "w", "sigma.p", "sigma.lam", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "alpha.p", "theta", "w", "sigma.p", "sigma.lam", "fit", "fit.new")])

stopCluster(cl)


#----------Overdispersion in both abund and detection + kriging + min covariates------------
# Define model
sink("Code/Jags_Models/sp_min_od_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
beta1.lam ~ dnorm(0, 0.01)
beta2.lam ~ dnorm(0, 0.01)
    
    
    alpha.p ~ dnorm(0, 0.01)
 beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[], Omega[,]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    # theta ~ dgamma(0.01, 0.01) # spatial correlation rate
    # theta ~ dunif(0.0001, 1000) # spatial correlation rate
    theta ~ dnorm(0, 0.01)T(0.00001, 1000) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
    mu.lam[i] ~ dnorm(alpha.lam, tau.lam)
    # muW[i] <- mu.lam[i]
    muW[i] <- 0
    }
    
    tau.lam <- 1 / (sigma.lam * sigma.lam)
    sigma.lam ~ dunif(0, 3)
    
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + mu.lam[i] + W[i] # + beta2.lam*elev[i]*elev[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       beta2.lam = rnorm(1, 0, 1),
       # beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       theta = runif(1, 2, 50),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam",
             "beta1.lam",
             "beta2.lam",
             "theta",
             "w",
             "alpha.p", 
             # "eps.lam",
             "sigma.lam",
             "beta10.p",
             "delta.p",
             # "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     # elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover * gcover,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/sp_min_od_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta10.p", "theta", "w", "sigma.p", "sigma.lam", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta10.p", "theta", "w", "sigma.p", "sigma.lam", "fit", "fit.new")])

stopCluster(cl)



#----------Overdispersion in both abund and detection + kriging + simple covariates----No--------
# Define model
sink("Code/Jags_Models/simple_sp_od_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    # beta4.lam ~ dnorm(0, 0.01)
    # beta5.lam ~ dnorm(0, 0.01)
    # beta6.lam ~ dnorm(0, 0.01)
    # beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    # beta10.lam ~ dnorm(0, 0.01)
    # beta11.lam ~ dnorm(0, 0.01)
    beta12.lam ~ dnorm(0, 0.01)
    # beta13.lam ~ dnorm(0, 0.01)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[], Omega[,]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    # theta ~ dgamma(0.01, 0.01) # spatial correlation rate
# theta ~ dunif(0.0001, 1000) # spatial correlation rate
    theta ~ dnorm(0, 0.01)T(0.00001, 1000) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
mu.lam[i] ~ dnorm(alpha.lam, tau.lam)
    muW[i] <- mu.lam[i]
    }

tau.lam <- 1 / (sigma.lam * sigma.lam)
sigma.lam ~ dunif(0, 3)

    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta9.lam*canopy[i] + W[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       # beta4.lam = rnorm(1, 0, 1),
       # beta5.lam = rnorm(1, 0, 1),
       # beta6.lam = rnorm(1, 0, 1),
       # beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       # beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       theta = runif(1, 2, 50),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

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
             "sigma.lam",
             "delta.p",
             # "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     # elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover * gcover,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/simple_sp_od_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "theta", "w", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)




#----------Overdispersion in detection + kriging + simple covariates------------
# Define model
sink("Code/Jags_Models/simple_sp_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    # beta4.lam ~ dnorm(0, 0.01)
    # beta5.lam ~ dnorm(0, 0.01)
    # beta6.lam ~ dnorm(0, 0.01)
    # beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    # beta10.lam ~ dnorm(0, 0.01)
    # beta11.lam ~ dnorm(0, 0.01)
    beta12.lam ~ dnorm(0, 0.01)
    # beta13.lam ~ dnorm(0, 0.01)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[], Omega[,]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    # theta ~ dgamma(0.01, 0.01) # spatial correlation rate
    theta ~ dnorm(0, 0.01)T(0.00001, 1000) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
  
    muW[i] <- 0
    }
    
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta9.lam*canopy[i] + W[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       # beta4.lam = rnorm(1, 0, 1),
       # beta5.lam = rnorm(1, 0, 1),
       # beta6.lam = rnorm(1, 0, 1),
       # beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       # beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       theta = runif(1, 2, 50),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

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
             # "sigma.lam",
             "delta.p",
             # "sigma.site",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     # elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover * gcover,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/simple_sp_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "theta", "w", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "fit", "fit.new")])

stopCluster(cl)




#----------Overdispersion in detection simple spatial model no elevation^2---YES? theta and W weird?----------
# Define model
sink("Code/Jags_Models/elev_sp_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    # beta4.lam ~ dnorm(0, 0.01)
    # beta5.lam ~ dnorm(0, 0.01)
    # beta6.lam ~ dnorm(0, 0.01)
    # beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    # beta10.lam ~ dnorm(0, 0.01)
    # beta11.lam ~ dnorm(0, 0.01)
    beta12.lam ~ dnorm(0, 0.01)
    # beta13.lam ~ dnorm(0, 0.01)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[1:n.transects], Omega[1:n.transects, 1:n.transects]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    # theta ~ dgamma(0.01, 0.01) # spatial correlation rate
# theta ~ dunif(0, 1000) # spatial correlation rate
theta ~ dnorm(0, 0.01) # T(0, 1000) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
    # muW[i] <- alpha.lam 
# muW[i] <- 0
    }
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
muW[1:n.transects] <- 0

    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    # log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta8.lam*ltwi[i] + beta9.lam*canopy[i] + W[i]

    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta8.lam*ltwi[i] + beta9.lam*canopy[i] + W[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       # beta4.lam = rnorm(1, 0, 1),
       # beta5.lam = rnorm(1, 0, 1),
       # beta6.lam = rnorm(1, 0, 1),
       # beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       # beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       theta = runif(1, 2, 50),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

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

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     # elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover * gcover,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d_km,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/elev_sp_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "theta", "w", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max







#----------Full model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/full_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max





#----------Full model only elev site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/full_elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)#T(-3, 5)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)

    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("sigma.site")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max


#----------Full model only elev site abund Overdispersion in detection *uniform priors*------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/full_elev_od_unif.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dunif(-5, 5)
    beta1.lam ~ dunif(-5, 5)
    beta2.lam ~ dunif(-5, 5)
    beta3.lam ~ dunif(-5, 5)
    beta4.lam ~ dunif(-5, 5)
    beta5.lam ~ dunif(-5, 5)
    beta6.lam ~ dunif(-5, 5)
    beta7.lam ~ dunif(-5, 5)
    beta8.lam ~ dunif(-5, 5)
    beta9.lam ~ dunif(-5, 5)
    beta10.lam ~ dunif(-5, 5)
    beta11.lam ~ dunif(-5, 5)
    beta12.lam ~ dunif(-5, 5)
    beta13.lam ~ dunif(-5, 5)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dunif(-5, 5)
    beta1.p ~ dunif(-5, 5)
    beta2.p ~ dunif(-5, 5)
    beta3.p ~ dunif(-5, 5)
    beta4.p ~ dunif(-5, 5)
    beta5.p ~ dunif(-5, 5)
    beta10.p ~ dunif(-5, 5)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = runif(1, -1, 1),
       beta1.lam = runif(1, -1, 1),
       beta2.lam = runif(1, -1, 1),
       beta3.lam = runif(1, -1, 1),
       beta4.lam = runif(1, -1, 1),
       beta5.lam = runif(1, -1, 1),
       beta6.lam = runif(1, -1, 1),
       beta7.lam = runif(1, -1, 1),
       beta8.lam = runif(1, -1, 1),
       beta9.lam = runif(1, -1, 1),
       beta10.lam = runif(1, -1, 1),
       beta11.lam = runif(1, -1, 1),
       beta12.lam = runif(1, -1, 1),
       beta13.lam = runif(1, -1, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = runif(1, -1, 1),
       beta1.p = runif(1, -1, 1),
       beta2.p = runif(1, -1, 1),
       beta3.p = runif(1, -1, 1),
       beta4.p = runif(1, -1, 1),
       beta5.p = runif(1, -1, 1),
       beta10.p = runif(1, -1, 1),
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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_elev_od_unif.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)





#----------min elev2 model site abund no od in detection----Doesn't Fit--------
# random effect of site on abundance in detection
# Define model
sink("Code/Jags_Models/min.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
       sigma.site = runif(1, 0, 1))#,
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
             # "delta.p",
             "sigma.site",
             # "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/min.txt", pjor.od.data, inits, n.adapt=3000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)




#----------min elev2 model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/min_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/min_od.txt", pjor.od.data, inits, n.adapt=3000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max







#----------min elev2 model Overdispersion in abundance and detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/min_od_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    eps.lam[i] ~ dnorm(0, tau.site)
    }

    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)

    for(i in 1:n.transects){
    trans.lam[i] ~ dnorm(alpha.lam, tau.trans)
    }

    sigma.trans ~ dunif(0, 3)
    tau.trans <- 1/(sigma.trans*sigma.trans)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + eps.lam[site[i]] + trans.lam[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
       sigma.trans = runif(1, 0, 1),
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
             "sigma.trans",
             "sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/min_od_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max











#----------min elev model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/min_elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/min_elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

par(mfrow = c(1,1))
plot((sd(Data$elev)*pjor.od.data$elev)+mean(Data$elev), jitter(pjor.od.data$C[ , 1]))
for(i in 2:6) {
  points((sd(Data$elev)*pjor.od.data$elev)+mean(Data$elev), jitter(pjor.od.data$C[ , i]))
}

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max













#----------min elev2 and stream model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/min_stream_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta13.lam*stream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
                     stream = stream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/min_stream_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max








#----------Full model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/full_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max





#----------Full model only elev site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/full_elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)#T(-3, 5)
    beta2.lam ~ dnorm(0, 0.01)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/full_elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("sigma.site")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max



#---------- extra? ---------------



#----------Overdispersion in abundance and detection linear elev -okay converge, decent fit -- Use this model---------------
# Define model
sink("pjor_od3.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    #beta2.lam ~ dnorm(0, 0.01)
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
    eps.lam[i] ~ dnorm(0, tau.site)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta2.lam*elev2[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
             #"beta2.lam", 
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
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     #elev2 = elev2, 
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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od3.txt", pjor.od.data, inits, n.adapt=nb00, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=nt00) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od3<- mcmc.list(out)
plot(pjor_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor3 <- mean(pjor_od3[, "fit.new",][[i]] > pjor_od3[, "fit",][[i]]) # 0.31 good 
print(bayesP.pjor3, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od3[, "fit",]), as.matrix(pjor_od3[, "fit.new",])) # imperfect but acceptable123
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # okay

N.pjor <- matrix(NA, 195, 3)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od3[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.pjor) <- c("CI_2.5", "Median", "CI_97.5")
N.pjor
cbind(NaiveRichness, PJORmin, N.pjor)

# find percent greater than zero
Psi.jor <- NA # probability that site is occupied
for(i in 1:195){
  Psi.jor[i] <- length(which(as.numeric(as.matrix(pjor_od3[, c(paste("N[", i, "]", sep = ""))])) > 0))/3000
}

print(Psi.jor, dig = 3)

print(cbind(NaiveRichness, PJORmin, Psi.jor), dig = 3)

#----------Overdispersion in abundance and detection no slope2 -no converge--------------
# Define model
sink("pjor_od4.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta10.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    #beta12.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

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
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, -.5, 1),
       beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       #beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, -.5, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -.5, 1),
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
             #"beta12.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
                     slope = slope,
                     #slope2 = slope2, 
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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(3)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od4.txt", pjor.od.data, inits, n.adapt=nb00, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=nt00) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od4<- mcmc.list(out)
plot(pjor_od4[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od4[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor3 <- mean(pjor_od3[, "fit.new",][[i]] > pjor_od3[, "fit",][[i]]) # 0.31 good 
print(bayesP.pjor3, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od3[, "fit",]), as.matrix(pjor_od3[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od4[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])


N.pjor <- matrix(NA, 195, 3)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od3[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.pjor) <- c("CI_2.5", "Median", "CI_97.5")
N.pjor
cbind(NaiveRichness, PJORmin, N.pjor)

# find percent greater than zero
Psi.jor <- NA # probability that site is occupied
for(i in 1:195){
  Psi.jor[i] <- length(which(as.numeric(as.matrix(pjor_od3[, c(paste("N[", i, "]", sep = ""))])) > 0))/3000
}

print(Psi.jor, dig = 3)

print(cbind(NaiveRichness, PJORmin, Psi.jor), dig = 3)


#----------Simple model with overdispersion in detection- okay not great------------
# Define model
sink("pjor_od4.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    #beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    #beta4.lam ~ dnorm(0, 0.01)
    #beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta10.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    #beta12.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       #beta4.lam = rnorm(1, 0, 1),
       #beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       #beta12.lam = rnorm(1, 0, 1),
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
             "beta3.lam",
             #"beta4.lam",
             #"beta5.lam",
             "beta6.lam",
             "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             "beta10.lam",
             "beta11.lam",
             #"beta12.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     #elev2 = elev2, 
                     slope = slope,
                     #slope2 = slope2, 
                     #aspectN = aspectN,
                     #aspectE = aspectE,
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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od4.txt", pjor.od.data, inits, n.adapt=300000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od4 <- mcmc.list(out)
plot(pjor_od4[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od3[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor4 <- mean(pjor_od4[, "fit.new",][[i]] > pjor_od4[, "fit",][[i]]) # 0.295
print(bayesP.pjor4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od4[, "fit",]), as.matrix(pjor_od4[, "fit.new",])) # okay
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_od4[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # 



#
#----------Overdispersion in abundance and detection linear elev -okay converge, decent fit -- Use this model---------------
# Define model
sink("pjor_od5.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
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
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
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
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
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
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od5.txt", pjor.od.data, inits, n.adapt=nb00, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=nt00) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od5<- mcmc.list(out)
plot(pjor_od5[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od5[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor5 <- mean(pjor_od5[, "fit.new",][[i]] > pjor_od3[, "fit",][[i]]) # 
print(bayesP.pjor5, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od5[, "fit",]), as.matrix(pjor_od5[, "fit.new",])) #
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_od5[,c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), dig=3) # okay

N.pjor <- matrix(NA, 195, 3)
for(i in 1:195){
  N.pjor[i, ] <- apply(as.matrix(pjor_od5[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.pjor) <- c("CI_2.5", "Median", "CI_97.5")
N.pjor
cbind(NaiveRichness, PJORmin, N.pjor)

# find percent greater than zero
Psi.jor <- NA # probability that site is occupied
for(i in 1:195){
  Psi.jor[i] <- length(which(as.numeric(as.matrix(pjor_od5[, c(paste("N[", i, "]", sep = ""))])) > 0))/3000
}

print(Psi.jor, dig = 3)

print(cbind(NaiveRichness, PJORmin, Psi.jor), dig = 3)






#----------Overdispersion in detection simple spatial model - elev^2 doesn't converge --- -------------

# Define model
sink("pjor_od5.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    # beta4.lam ~ dnorm(0, 0.01)
    # beta5.lam ~ dnorm(0, 0.01)
    # beta6.lam ~ dnorm(0, 0.01)
    # beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    # beta10.lam ~ dnorm(0, 0.01)
    # beta11.lam ~ dnorm(0, 0.01)
    beta12.lam ~ dnorm(0, 0.01)
    # beta13.lam ~ dnorm(0, 0.01)
    
    # for(i in 1:n.sites){
    # eps.lam[i] ~ dnorm(0, tau.site)
    # }
    # 
    # sigma.site ~ dunif(0, 3)
    # tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 3)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Spatial component (kriging)
    W[1:n.transects] ~ dmnorm(muW[], Omega[,]) 
    tauw ~ dgamma(0.001, 0.001)
    w <- 1/sqrt(tauw) # spatial SD
    theta ~ dgamma(0.01, 0.01) # spatial correlation rate
    for (i in 1:n.transects){ 
    for(j in 1:n.transects){
    # H[i, j] <- (1/tauw) * exp(-theta*pow(d[i,j] ,2)) # Spatial covariance matrix
    # H[i, j] <- (1/tauw) * exp(-theta*d[i,j]) # Spatial covariance matrix
    H[i, j] <- (1/tauw) * exp(-d[i,j] / theta) # Spatial covariance matrix
    }
    muW[i] <- 0
    }
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta12.lam*slope2[i] + W[i] # + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

Nst <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 2, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       # beta4.lam = rnorm(1, 0, 1),
       # beta5.lam = rnorm(1, 0, 1),
       # beta6.lam = rnorm(1, 0, 1),
       # beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       # beta13.lam = rnorm(1, 0, 1),
       # eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       # sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
             "beta3.lam",
             # "beta4.lam",
             # "beta5.lam",
             # "beta6.lam",
             # "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             # "beta10.lam",
             # "beta11.lam",
             "beta12.lam",
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

pjor.od.data <- list(C = as.matrix(PJOR[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
                     slope = slope,
                     slope2 = slope2, 
                     # aspectN = aspectN,
                     # aspectE = aspectE,
                     ltwi = ltwi,
                     # tpi = tpi,
                     trail = trail,
                     canopy = canopy,
                     gcover = gcover,
                     gcover2 = gcover2,
                     litterdepth = litterdepth,
                     # lstream = lstream,
                     site = as.numeric(site),
                     Temp.s = as.matrix(Temp.s[ ,1:6]),
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     d = d,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od5.txt", pjor.od.data, inits, n.adapt=100, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od2.txt", pjor.od.data, inits, n.adapt=nb00, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=nt00) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od<- mcmc.list(out)
plot(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta8.lam", "beta9.lam", "beta12.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "theta", "w", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

plot(pjor_od[[1]][,c("fit.new")]) # problem chain
plot(pjor_od[[2]][,c("fit.new")])
plot(pjor_od[[3]][,c("fit.new")])
plot(pjor_od[[4]][,c("fit.new")]) # 

pjor_od <- mcmc.list(out[[2]], out[[3]], out[[3]])
# Check fit
for(i in 1:3) bayesP.pjor <- mean(pjor_od[, "fit.new",][[i]] > pjor_od[, "fit",][[i]]) # 0.332 good fit but problems with convergence of elev and elev^2
print(bayesP.pjor, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_od[, "fit",]), as.matrix(pjor_od[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])
print(gelman.diag(x=pjor_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # no convergence of elev and elev^2
g <- matrix(NA, nrow=nvar(pjor_od), ncol=2)
for (v in 1:nvar(pjor_od)) {
  g[v,] <- gelman.diag(pjor_od[,v])$psrf
}
max(g[,1], na.rm = TRUE) # a bit higher than multivariate 
max(g[,2], na.rm = TRUE) # max





#----------PJOR Occupancy Model----rerun------
# Define model
sink("pjor_od1_occ.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    #beta2.lam ~ dnorm(0, 0.01)
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
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    sigma.site ~ dunif(0, 25)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 25)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    z[i] ~ dbern(psi[i])
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    lpsi.lim[i] <- min(999, max(-999, lpsi[i])) # Help stabilize the logit
    psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))    
    
    for(j in 1:n.surveys){
    y[i, j] ~ dbern(eff.p[i, j])
    eff.p[i,j] <- z[i] * p[i,j]
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    Presi[i,j] <- abs(y[i,j] - p[i,j]) # Absolute residual
    y.new[i,j] ~ dbern(eff.p[i,j])
    Presi.new[i,j] <- abs(y.new[i,j] - p[i,j])
    } 
    }
    fit <- sum(Presi[,]) # Discrepancy for actual data set
    fit.new <- sum(Presi.new[,]) #Discrepancy for replicate data set
    }
    ", fill = TRUE)
sink()

# Convert EWIL Data to 0/1
PJOR01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    PJOR01[i,j] <- ifelse(PJOR[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(PJOR01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       #beta2.lam = rnorm(1, 0, 1),
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
             #"beta2.lam", 
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
             #"eps.lam",
             # "delta.p",
             "sigma.site",
             "sigma.p",
             "psi",
             "z",
             #"p",
             "fit",
             "fit.new")

pjor.od.data01 <- list(y = as.matrix(PJOR01[, 1:6]), 
                       n.transects = n.transects, 
                       n.surveys = n.surveys,
                       n.sites = length(site.inits), 
                       elev = elev, 
                       #elev2 = elev2, 
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
                       Temp.s = Temp.s[ ,1:6],
                       Temp.s2 = Temp.s*Temp.s,
                       RH.s = RH.s[ ,1:6],
                       Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od1_occ.txt", pjor.od.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_occ_od <- mcmc.list(out)
plot(pjor_occ_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_occ_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor_occ_od <- mean(pjor_occ_od[, "fit.new",][[i]] > pjor_occ_od[, "fit",][[i]]) # 0.442 good 
print(bayesP.pjor_occ_od, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_occ_od[, "fit",]), as.matrix(pjor_occ_od[, "fit.new",])) # imperfect but acceptable123
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_occ_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), dig=3) # 


#----------PJOR Occupancy Model----doesn't work with site effect------
# Define model
sink("pjor_occ.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    #beta2.lam ~ dnorm(0, 0.01)
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
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    sigma.site ~ dunif(0, 15)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    # Likelihood
    for(i in 1:n.transects){
    z[i] ~ dbern(psi[i])
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    lpsi.lim[i] <- min(999, max(-999, lpsi[i])) # Help stabilize the logit
    psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))    
    
    for(j in 1:n.surveys){
    y[i, j] ~ dbern(eff.p[i, j])
    eff.p[i,j] <- z[i] * p[i,j]
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    Presi[i,j] <- abs(y[i,j] - p[i,j]) # Absolute residual
    y.new[i,j] ~ dbern(eff.p[i,j])
    Presi.new[i,j] <- abs(y.new[i,j] - p[i,j])
    } 
    }
    fit <- sum(Presi[,]) # Discrepancy for actual data set
    fit.new <- sum(Presi.new[,]) #Discrepancy for replicate data set
    }
    ", fill = TRUE)
sink()

# Convert EWIL Data to 0/1
PJOR01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    PJOR01[i,j] <- ifelse(PJOR[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(PJOR01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       #beta2.lam = rnorm(1, 0, 1),
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
       #eps.lam = rnorm(48, 0, 10),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 5, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             #"beta2.lam", 
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
             #"eps.lam",
             # "delta.p",
             "sigma.site",
             #"sigma.p",
             "psi",
             "z",
             #"p",
             "fit",
             "fit.new")

pjor.data01 <- list(y = as.matrix(PJOR01[, 1:6]), 
                    n.transects = n.transects, 
                    n.surveys = n.surveys,
                    n.sites = length(site.inits), 
                    elev = elev, 
                    #elev2 = elev2, 
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
                    Temp.s = Temp.s[ ,1:6],
                    Temp.s2 = Temp.s*Temp.s,
                    RH.s = RH.s[ ,1:6],
                    Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_occ.txt", pjor.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni00, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_occ <- mcmc.list(out)
plot(pjor_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.pjor_occ <- mean(pjor_occ_od[, "fit.new",][[i]] > pjor_occ_od[, "fit",][[i]]) #
print(bayesP.pjor_occ_od, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(pjor_occ[, "fit",]), as.matrix(pjor_occ[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=pjor_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # 




#----------Stream x Elevation^2 model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/stream-elev2_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 5)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta6.lam*stream[i] + beta8.lam*stream[i]*elev[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()


#----------Stream x Elevation model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/stream-elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 5)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta6.lam*stream[i] + beta8.lam*streamelev[i] + beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

#----------Final model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/final_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, 1)
    }
    
    sigma.site ~ dt(0, 1 / (25^2), 1)I(0, ) 	## implies half-cauchy with scale of 25 which is probably absurdly large
    # tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
# Z.lam ~ dnorm(0, 1)
# Z.p ~ dnorm(0, 1)

    sigma.p ~ dt(0, 1 / (25^2), 1)I(0, ) 	## implies half-cauchy with scale of 25
    # tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, 1)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*twi[i] + beta4.lam*litterdepth[i] + beta5.lam*gcover[i] + beta6.lam*stream[i] + sigma.site * eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + sigma.p * delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

#----------Final model site abund Overdispersion in detection uniform------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
# Define model
sink("Code/Jags_Models/final_od_unif.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dunif(-12, 12)
    beta1.lam ~ dunif(-12, 12)
    beta2.lam ~ dunif(-12, 12)
    beta3.lam ~ dunif(-12, 12)
    beta4.lam ~ dunif(-12, 12)
    beta5.lam ~ dunif(-12, 12)
    beta6.lam ~ dunif(-12, 12)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dunif(-12, 12)
    beta1.p ~ dunif(-12, 12)
    beta2.p ~ dunif(-12, 12)
    beta3.p ~ dunif(-12, 12)
    beta4.p ~ dunif(-12, 12)
    beta5.p ~ dunif(-12, 12)
    beta10.p ~ dunif(-12, 12)
    
    sigma.p ~ dunif(0, 10)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*twi[i] + beta4.lam*litterdepth[i] + beta5.lam*gcover[i] + beta6.lam*stream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()



#----------Final model site abund Overdispersion in detection------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/final_elev_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    # beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 5)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*twi[i] + beta4.lam*litterdepth[i] + beta5.lam*gcover[i] + beta6.lam*stream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()

#----------Final model site abund Overdispersion in detection elev^2 and stream^2------------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("Code/Jags_Models/final_stream2_od.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 5)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(alpha.p, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*twi[i] + beta4.lam*litterdepth[i] + beta5.lam*gcover[i] + beta6.lam*stream[i] + beta7.lam*stream[i]*stream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <-  beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    }
    }
    
    # Derived quantities
    totalN<- sum(N[1:n.transects])
    fit <- sum(E[1:n.transects, 1:n.surveys])
    fit.new <- sum(E.new[1:n.transects, 1:n.surveys])
    
    }
    ", fill = TRUE)
sink()