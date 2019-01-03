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



# check sorting and organization of elevation data

# just use temperature or precip in place of elevation

str(Data)

str(df)




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
                     tmean = df$tmean,
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
