#########################################################################
# GSMNP Elevation Gradient Analysis
# October 31, 2012
# Daniel J. Hocking
########################################################################

#-------Set WD and Load Packages------
# setwd('/Users/Dan/Documents/Research/SmokysElevation/Analysis/Bayesian/')

library(parallel)
library(rjags)
library(dplyr)
library(ggplot2)

if(!dir.exists("Results/JAGS")) dir.create("Results/JAGS", recursive = TRUE)

######

#-----------Load Data and do Manipulation--------
# Import and check data
Data <- read.table('Data/Raw/TransectData2012.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)

Side <- read.table('Data/Raw/Side.csv', header = TRUE, sep = ',')
Data <- merge(x=Data, y=Side, by.x = "site", by.y = 'site', all.x = TRUE, sort = FALSE)

Stream <- read.table('Data/Raw/Strm_Dist.csv', header = TRUE, sep = ",")
Data <- merge(x = Data, y=Stream, by.x = "transect", by.y = 'plot_trans', all.x = TRUE, sort = FALSE)

Data <- Data[order(Data$elev, na.last = TRUE), ] 


PJOR <- Data[ , 3:8]
Precip <- Data[ , c('precip1', 'precip2', 'precip3', 'precip4', 'precip5', 'precip6')]
Temp <- Data[ , c('temp1', 'temp2', 'temp3', 'temp4', 'temp5', 'temp6')]
RH <- Data[ , c('RH1', 'RH2', 'RH3', 'RH4', 'RH5', 'RH6')]
Rain <- Data[ , c('rain1', 'rain2', 'rain3', 'rain4', 'rain5', 'rain6')]
VegWet <- Data[ , c('veg1', 'veg2', 'veg3', 'veg4', 'veg5', 'veg6')]
Ground <- Data[ , c('ground1', 'ground2', 'ground3', 'ground4', 'ground5', 'ground6')]
DCON <- Data[ , c('DCON1', 'DCON2', 'DCON3', 'DCON4', 'DCON5', 'DCON6')]
DIMI <- Data[ , c('DIMI1', 'DIMI2', 'DIMI3', 'DIMI4', 'DIMI5', 'DIMI6')]
DMON <- Data[ , c('DMON1', 'DMON2', 'DMON3', 'DMON4', 'DMON5', 'DMON6')]
DOCO <- Data[ , c('DOCO1', 'DOCO2', 'DOCO3', 'DOCO4', 'DOCO5', 'DOCO6')]
DQUA <- Data[ , c('DQUA1', 'DQUA2', 'DQUA3', 'DQUA4', 'DQUA5', 'DQUA6')]
DSAN <- Data[ , c('DSAN1', 'DSAN2', 'DSAN3', 'DSAN4', 'DSAN5', 'DSAN6')]
DWRI <- Data[ , c('DWRI1', 'DWRI2', 'DWRI3', 'DWRI4', 'DWRI5', 'DWRI6')]
EWIL <- Data[ , c('EWIL1', 'EWIL2', 'EWIL3', 'EWIL4', 'EWIL5', 'EWIL6')]
GPOR <- Data[ , c('GPOR1', 'GPOR2', 'GPOR3', 'GPOR4', 'GPOR5', 'GPOR6')]
NVIR <- Data[ , c('NVIR1', 'NVIR2', 'NVIR3', 'NVIR4', 'NVIR5', 'NVIR6')]
PGLU <- Data[ , c('PGLU1', 'PGLU2', 'PGLU3', 'PGLU4', 'PGLU5', 'PGLU6')]
PSER <- Data[ , c('PSER1', 'PSER2', 'PSER3', 'PSER4', 'PSER5', 'PSER6')]
PTEY <- Data[ , c('PTEY1', 'PTEY2', 'PTEY3', 'PTEY4', 'PTEY5', 'PTEY6')]

DCON <- replace(DCON, is.na(PJOR), NA)
DIMI <- replace(DIMI, is.na(PJOR), NA)
DMON <- replace(DMON, is.na(PJOR), NA)
DOCO <- replace(DOCO, is.na(PJOR), NA)
DQUA <- replace(DQUA, is.na(PJOR), NA)
DSAN <- replace(DSAN, is.na(PJOR), NA)
DWRI <- replace(DWRI, is.na(PJOR), NA)
EWIL <- replace(EWIL, is.na(PJOR), NA)
GPOR <- replace(GPOR, is.na(PJOR), NA)
NVIR <- replace(NVIR, is.na(PJOR), NA)
PGLU <- replace(PGLU, is.na(PJOR), NA)
PSER <- replace(PSER, is.na(PJOR), NA)
PTEY <- replace(PTEY, is.na(PJOR), NA)

# Total Captures
Count <- c( 
sum(DCON, na.rm = TRUE),
sum(DIMI, na.rm = TRUE),
sum(DMON, na.rm = TRUE),
sum(DOCO, na.rm = TRUE),
sum(DQUA, na.rm = TRUE),
sum(DSAN, na.rm = TRUE),
sum(DWRI, na.rm = TRUE),
sum(EWIL, na.rm = TRUE),
sum(GPOR, na.rm = TRUE),
sum(NVIR, na.rm = TRUE),
sum(PGLU, na.rm = TRUE),
sum(PJOR, na.rm = TRUE),
sum(PSER, na.rm = TRUE),
sum(PTEY, na.rm = TRUE))

Species <- c('D. conanti', 'D. imitator', 'D. monticola', 'D. ocoee', 'D. quadramaculatus', 'D. santeetlah', 'D. wrighti', 'E. wilderae', 'G. porphyriticus', 'N. viridescens', 'P. glutinosus', 'P. jordani', 'P. seratus', 'P. teyahalee')

write.table(data.frame(Species, Count), file = 'Data/TotalCounts.csv', sep = ",", row.names = FALSE)

PJORmin <- apply(PJOR, 1, function(x) max(x, na.rm = TRUE))
DCONmin <- apply(DCON, 1, function(x) max(x, na.rm = TRUE))
DIMImin <- apply(DIMI, 1, function(x) max(x, na.rm = TRUE))
DMONmin <- apply(DMON, 1, function(x) max(x, na.rm = TRUE))
DOCOmin <- apply(DOCO, 1, function(x) max(x, na.rm = TRUE))
DQUAmin <- apply(DQUA, 1, function(x) max(x, na.rm = TRUE))
DSANmin <- apply(DSAN, 1, function(x) max(x, na.rm = TRUE))
DWRImin <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE))
EWILmin <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE))
GPORmin <- apply(GPOR, 1, function(x) max(x, na.rm = TRUE))
NVIRmin <- apply(NVIR, 1, function(x) max(x, na.rm = TRUE))
PGLUmin <- apply(PGLU, 1, function(x) max(x, na.rm = TRUE))
PSERmin <- apply(PSER, 1, function(x) max(x, na.rm = TRUE))
PTEYmin <- apply(PTEY, 1, function(x) max(x, na.rm = TRUE))

Naive <- data.frame(cbind(DCONmin, DIMImin, DMONmin, DOCOmin, DQUAmin, DSANmin, DWRImin, EWILmin, GPORmin, NVIRmin, PGLUmin, PJORmin, PSERmin, PTEYmin))

summary(Naive)

NaiveRichness <- data.frame(1:195, Data$transect, Data[,"side"], Data[,"elev"], apply(as.matrix(Naive), 1, function(x) sum(x > 0)))
colnames(NaiveRichness) <- c('transect.num', 'transect', 'side', 'elev', 'richness')

plot(NaiveRichness$elev, NaiveRichness$richness)
lines(smooth.spline(NaiveRichness$elev, NaiveRichness$richness, all.knots = FALSE, nknots = 20))
lines(lowess(NaiveRichness$elev, NaiveRichness$richness), col = 'red')

lmPJOR <- lm(PJORmin ~ NaiveRichness$elev + I(NaiveRichness$elev^2))
lmPTEY <- lm(PTEYmin ~ NaiveRichness$elev + I(NaiveRichness$elev^2))

plot(NaiveRichness$elev, PJORmin, type = "p", xlab = 'Elevation (m)', ylab = 'Max count of Plethodon sp.')
points(NaiveRichness$elev, PTEYmin, col = 2)
#lines(lowess(NaiveRichness$elev, PJORmin), col = 1)
#lines(smooth.spline(NaiveRichness$elev, PJORmin, nknots = 10), col = 'red')
lines(Data$elev, fitted(lmPJOR), col = 1)
lines(Data$elev, fitted(lmPTEY), col = 2)

lmDWRI <- lm(DWRImin ~ NaiveRichness$elev + I(NaiveRichness$elev^2))
lmDOCO <- lm(DOCOmin ~ NaiveRichness$elev + I(NaiveRichness$elev^2))

plot(NaiveRichness$elev, DWRImin, type = "p", xlab = 'Elevation (m)', ylab = 'Max count of Desmognathus')
points(NaiveRichness$elev, DOCOmin, col = 2)
#lines(lowess(NaiveRichness$elev, PJORmin), col = 1)
#lines(smooth.spline(NaiveRichness$elev, PJORmin, nknots = 10), col = 'red')
lines(Data$elev, fitted(lmDWRI), col = 1)
lines(Data$elev, fitted(lmDOCO), col = 2)

n.transects <- length(Data$gcover)

Data <- Data %>%
  dplyr::arrange(site)

Data$site1 <- as.numeric(seq(1,n.transects))
for(i in 2:n.transects){
  Data$site1[i] <- ifelse(Data$site[i] == Data$site[i-1], Data$site1[i-1], (Data$site1[i-1] + 1))
}
Data$sitef <- as.factor(Data$site1)

Data$sidef <- as.factor(Data$side)

# Prepare and standardize data
precip.m <- mean(as.matrix(Precip), na.rm = TRUE)
precip.sd <- sd(unlist(as.vector(Precip)), na.rm = TRUE)

Precip.s <- as.matrix((Precip - precip.m)/precip.sd) # Scale over all precip
mean(as.matrix(Precip.s), na.rm = TRUE) # check that = 0
sd(unlist(as.vector(Precip.s)), na.rm = TRUE) # check that = 1
dimnames(Precip.s) <- NULL
str(Precip.s)
dim(Precip.s)

temp.m <- mean(as.matrix(Temp), na.rm = TRUE)
temp.sd <- sd(unlist(as.vector(Temp)), na.rm = TRUE)
Temp.s <- as.matrix((Temp - temp.m)/temp.sd) # Scale over entire matrix
mean(as.matrix(Temp.s), na.rm = TRUE) # check that = 0
sd(unlist(as.vector(Temp.s)), na.rm = TRUE) # check that = 1
dimnames(Temp.s) <- NULL

RH.m <- mean(as.matrix(RH), na.rm = TRUE)
RH.sd <- sd(unlist(as.vector(RH)), na.rm = TRUE)
RH.s <- as.matrix((RH - RH.m)/RH.sd) # Scale over entire matrix
mean(as.matrix(RH.s), na.rm = TRUE) # check that = 0
sd(unlist(as.vector(RH.s)), na.rm = TRUE) # check that = 1
dimnames(RH.s) <- NULL

Temp2 <- Temp^2
Temp2.s <- (Temp2 - mean(Temp2, na.rm = TRUE))/sd(Temp2, na.rm = TRUE)

# PLJO <- as.matrix(PLJO)
# dimnames(PLJO) <- NULL

Rain <- as.matrix(Rain + 1)
dimnames(Rain) <- NULL

VegWet <- as.matrix(VegWet + 1)
dimnames(VegWet) <- NULL

Ground <- as.matrix(Ground + 1)
dimnames(Ground) <- NULL

tpi <- scale(as.vector(Data$tpi_100))
attributes(tpi) <- NULL
str(tpi)

solar <- (Data$solar - mean(Data$solar))/sd(Data$solar)
twi <- (Data$twi_10 - mean(Data$twi_10))/sd(Data$twi_10)
elev <- (Data$elev - mean(Data$elev))/sd(Data$elev)
Data$elev2 <- Data$elev^2
elev2 <- (Data$elev2 - mean(Data$elev2))/sd(Data$elev2)
slope <- (Data$Slope- mean(Data$Slope))/sd(Data$Slope)
Data$Slope2 <- Data$Slope^2
slope2 <- (Data$Slope2 - mean(Data$Slope2))/sd(Data$Slope2)
Data$AspectR <- Data$Aspect*2*pi/360 # convert degrees to radians
Data$AspectR2 <- Data$AspectR^2
Data$aspectN <- cos(Data$Aspect)
Data$aspectE <- sin(Data$Aspect)
aspectN <- Data$aspectN
aspectE <- Data$aspectE
trail <- Data$trail
site <- Data$sitef
canopy <- (Data$canopy - mean(Data$canopy))/sd(Data$canopy)
litterdepth <- (Data$litdepth - mean(Data$litdepth))/sd(Data$litdepth)
gcover <- (Data$gcover - mean(Data$gcover, na.rm = TRUE))/sd(Data$gcover, na.rm = TRUE)
# Data$ltwi <- log(Data$twi_10 + 10)
ltwi <- Data$twi_10
ltwi <- (Data$twi_10 - mean(Data$twi_10))/sd(Data$twi_10)
Data$lstream <- log(Data$strm_dist + 1)
lstream <- (Data$lstream - mean(Data$lstream))/sd(Data$lstream)
stream <- (Data$strm_dist - mean(Data$strm_dist))/sd(Data$ststrm_distream)
Data$gcover2 <- Data$gcover^2
gcover2 <- gcover^2

site.inits <- read.table('Data/Raw/site_inits.csv', header = TRUE, sep = ',')
site.inits <- as.vector(site.inits$siteinit)

n.transects <- nrow(Data)
n.surveys <- 6

# Scatterplot matrix
source('/Users/Dan/Documents/Statistics/R/Functions/panelcor.R')
Pairs1 <- data.frame(DWRImin, PJORmin, Data$elev, Data$tpi_100, Data$trail, Data$aspectN, Data$aspectE, Data$Slope, Data$solar, Data$twi_10)

pairs(Pairs1, upper.panel=panel.smooth, lower.panel=panel.cor) # nothing excessively correlated

Pairs2 <- data.frame(PJORmin, temp_30yr_active$tmean, ppt_30yr_active$ppt, Data$elev, Data$tpi_100, Data$Slope, Data$solar, Data$twi_10, canopy)

pairs(Pairs2, upper.panel=panel.smooth, lower.panel=panel.cor) # something not right!!!!!

pairs(select(df, -transect, -site, -Lat, -Long, -elev_bin, -lon, -lat, -trail), upper.panel=panel.smooth, lower.panel=panel.cor)

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

gcover <- impute.mean(gcover)
gcover2 <- gcover^2

gcover2 <- gcover^2
Temp.s2 <- Temp.s^2

VDamp <- matrix(0, 196, 6)
VWet <- matrix(0, 196, 6)
VDamp[VegWet == 2] <- 1
VWet[VegWet == 3] <- 1

Mist <- matrix(0, 196, 6)
Shower <- matrix(0, 196, 6)
Mist[Rain == 2] <- 1
Shower[Rain == 3] <- 1

# Summary
summary(Data)

Data$transectNum <- as.numeric(gsub(pattern="([0-9A-Za-z]+-)([0-9])", replacement="\\2", x=Data$transect, perl = FALSE))

range(aggregate(transectNum ~ site, data = Data, FUN = max)$transectNum) # number of transects per site

mean(aggregate(transectNum ~ site, data = Data, FUN = max)$transectNum) # number of transects per site

sum(aggregate(transectNum ~ site, data = Data, FUN = max)$transectNum) # total of 196 but only 195 used in analysis? Must be missing a transect (numbered but never run at 1 site)

UniqueSiteTransects <- aggregate(transectNum ~ site, data = Data, FUN = unique)
str(UniqueSiteTransects)
head(UniqueSiteTransects)
NumTrans <- 1:48
for(i in 1:length(UniqueSiteTransects$transectNum)){
  NumTrans[i] <- length(as.vector(UniqueSiteTransects[[2]][i][[1]]))
}
NumTrans
sum(NumTrans)

AgMax <- aggregate(transectNum ~ site, data = Data, FUN = max)
Test <- AgMax$transectNum == NumTrans

cbind(AgMax, NumTrans, Test, UniqueSiteTransects) # problem was no transect 4 at site 21B just 1,2,3 and 5,6,7

mean(NumTrans)

length(Data$transect) # 195

# find min and max elevation of each species
minMaxElev <- function(spp, elev = Data$elev) {
  elev.present <- NULL
  for(i in 1:dim(spp)[1]) {
    for(j in 1:dim(spp)[2]) {
      if(spp[i, j] > 0 & !is.na(spp[i, j])) {
        foo <- elev[i]
        elev.present <- c(elev.present, foo)
      }
    }
  }
  min.elev <- min(elev.present)
  max.elev <- max(elev.present)
  min.max.elev <- data.frame(min.elev, max.elev)
  return(min.max.elev)
}

minMaxElev(PJOR)

spp.list <- list(DCON=DCON, DIMI=DIMI, DOCO=DOCO, DMON=DMON, DQUA=DQUA, DSAN=DSAN, DWRI=DWRI, EWIL=EWIL, GPOR=GPOR, NVIR=NVIR, PGLU=PGLU, PJOR=PJOR, PSER=PSER, PTEY=PTEY)
min.max.elev <- t(sapply(spp.list, FUN = minMaxElev))
write.table(min.max.elev, file = "Data/min-max-elev.csv", sep = ",")


#---------------------------BAYESIAN----------------------------------

# remove NA in independent variables

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)) # should be a problem since no data in these dependent variables either
Temp.s <- impute.mean(Temp.s)
RH.s <- impute.mean(RH.s)
Precip.s <- impute.mean(Precip.s)

gcover <- impute.mean(gcover)

VegWet <- replace(VegWet, is.na(VegWet), 1)
Ground <- replace(Ground, is.na(Ground), 1)
Rain <- replace(Rain, is.na(Rain), 1)

VDamp <- matrix(0, 196, 6)
VWet <- matrix(0, 196, 6)
VDamp[VegWet == 2] <- 1
VWet[VegWet == 3] <- 1

Mist <- matrix(0, 196, 6)
Shower <- matrix(0, 196, 6)
Mist[Rain == 2] <- 1
Shower[Rain == 3] <- 1

#------------- PJOR max counts vs elev-------------

foo <- data.frame(site = Data$site, elevation = Data$elev, PJORmin)

ggplot(foo, aes(elevation, PJORmin)) + geom_point(aes(colour = site), position = "jitter") + geom_smooth()


#---------- set spatial matrix ---------

# set spatial relationships
library(sp)
library(rgdal)

# Data$transectf <- as.factor(Data$transect)

xy <- Data[ , c("Long", "Lat")] # add transect number?
names(xy) <- c("lon", "lat")
coordinates(xy) <- c("lon", "lat")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example

res <- spTransform(xy, CRS("+proj=utm +zone=17 ellps=NAD83"))
head(res)
str(res)

# plot_coords <- data.frame(transect = dimnames(res@coords)[[1]], res@coords, stringsAsFactors = FALSE)
# str(plot_coords)

d <- as.matrix(dist(sp::coordinates(res@coords[ , 1:2]))) # columns 1:2 are the x,y coordinates (UTM) - only problem with this method is keeping transect names. Will work if order doesn't change anywhere.

d <- d / 1000 # convert distances to kilometers

save.image("Data/jags_prep.RData")

################################################
# D. wrighti - Bayesian Analysis
################################################

#---------------Poisson D. Wrighti-*doesn't fit---- -----------------
# Define model
sink("dwri_p1.txt")
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j]
    
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

site.inits <- read.table('site_inits.csv', header = TRUE, sep = ',')
site.inits <- as.vector(site.inits$siteinit)

n.transects <- length(Data$elev)
n.surveys <- 6

Dwri.Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Dwri.Nst,
       alpha.lam = rnorm(1, 3, 1),
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
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1))#,
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
             "N",
             "p",
             "fit",
             "fit.new")

dwri.data <- list(C = as.matrix(DWRI[ , 1:6]), 
                  n.transects = length(elev), 
                  n.surveys = 6,
                  #n.sites = length(site.inits), 
                  elev = elev, 
                  elev2 = elev2, 
                  slope = slope,
                  slope2 = slope^2, 
                  aspectN = aspectN,
                  aspectE = aspectE,
                  ltwi = ltwi,
                  tpi = tpi,
                  trail = trail,
                  canopy = canopy,
                  gcover = gcover,
                  gcover2 = gcover^2,
                  litterdepth = litterdepth,
                  lstream = lstream,
                 # site = as.numeric(site),
                  Temp.s = as.matrix(Temp.s[ ,1:6]),
                  Temp.s2 = Temp.s*Temp.s,
                  RH.s = RH.s[ ,1:6],
                  Precip.s = Precip.s[ ,1:6])

saved.per.chain <- 1000
nc <- 3
nb <- 300000
nt <- 50
ni <- saved.per.chain*nt + nb

runif(1)

dwri.p1 <- jags(data = dwri.data, parameters.to.save = params, inits = inits, model.file = "dwri_p1.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)

# Posterior Predictive Checks
mean(dwri.p1$BUGSoutput$sims.list$fit.new > dwri.p1$BUGSoutput$sims.list$fit, na.rm = TRUE)
mean(dwri.p1$BUGSoutput$mean$fit) / mean(dwri.p1$BUGSoutput$mean$fit.new)
plot(dwri.p1$BUGSoutput$sims.list$fit, dwri.p1$BUGSoutput$sims.list$fit.new) # model doesn't fit
abline(0, 1, col = 'red')
**plot(dwri.p1$BUGSoutput$sims.list$alpha.lam, dwri.p1$BUGSoutput$sims.list$alpha.p)

traceplot(dwri.p1) # good mixing
par(mfrow = c(1,1))

# Evaluate convergence and if reasonable estimates
print(dwri.p1$BUGSoutput$summary[1:50, c(1,2,3, 4, 7,8,9)], dig = 3)
print(dwri.p1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(dwri.p1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(dwri.p1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(dwri.p1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)

(dwri.p1.tab <- data.frame(dwri.p1$BUGSoutput$summary[c("alpha.lam",
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
                                                         #"sigma.site",
                                                         "alpha.p",
                                                         "beta1.p",
                                                         "beta2.p",
                                                         "beta3.p",
                                                         "beta4.p",
                                                         "beta5.p",
                                                       # "beta6.p",
                                                        #"beta7.p",
                                                        #"beta8.p",
                                                        #"beta9.p",
                                                         "beta10.p"),
                                                       c(1,2,3,5,7,8,9)]))

# Check for autocorrelation in MCMC saved values
acf(dwri.p1$BUGSoutput$sims.list$alpha.lam)
acf(dwri.p1$BUGSoutput$sims.list$beta1.lam)
acf(dwri.p1$BUGSoutput$sims.list$beta2.lam)
acf(dwri.p1$BUGSoutput$sims.list$alpha.p)
acf(dwri.p1$BUGSoutput$sims.list$beta1.p)
acf(dwri.p1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(dwri.p1$BUGSoutput$sims.list$N[ ,1])
acf(dwri.p1$BUGSoutput$sims.list$N[ ,50])
acf(dwri.p1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(dwri.p1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(dwri.p1$BUGSoutput$sims.list$delta.p[ ,10,6])

#---------------Poisson D. Wrighti-FULL----poor fit moderate converge -----------------
sink("dwri_p2.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.001)
    beta1.lam ~ dnorm(0, 0.001)
    beta2.lam ~ dnorm(0, 0.001)
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
    
    alpha.p ~ dnorm(0, 0.001)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta6.p ~ dnorm(0, 0.01)
    beta7.p ~ dnorm(0, 0.01)
    beta8.p ~ dnorm(0, 0.01)
    beta9.p ~ dnorm(0, 0.001)
    beta10.p ~ dnorm(0, 0.01)
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + beta10.lam*gcover[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta6.p*VegDamp[i,j] + beta7.p*VegWet[i,j] + beta8.p*Mist[i,j] + beta9.p*Rain[i,j] + beta10.p*RH.s[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -2.4, 1),
       beta3.lam = rnorm(1, 0.9, 1),
       beta4.lam = rnorm(1, 0.2, .1),
       beta5.lam = rnorm(1, 0.3, .1),
       beta6.lam = rnorm(1, 0.5, .1),
       beta7.lam = rnorm(1, -0.5, 1),
       beta8.lam = rnorm(1, -0.1, .1),
       beta9.lam = rnorm(1, -0.15, .1),
       beta10.lam = rnorm(1, 0.1, 1),
       beta11.lam = rnorm(1, 0.2, 1),
       beta12.lam = rnorm(1, -1.1, 1),
       beta13.lam = rnorm(1, -.2, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, -0.2, .1),
       beta2.p = rnorm(1, -0.5, 1),
       beta3.p = rnorm(1, 0.2, .1),
       beta4.p = rnorm(1, 1, 1),
       beta5.p = rnorm(1, -0.5, 1),
       beta6.p = rnorm(1, 0.8, 1),
       beta7.p = rnorm(1, 1.5, 1),
       beta8.p = rnorm(1, 0, 1),
       beta9.p = rnorm(1, -0.5, 0.5),
       beta10.p = rnorm(1, 0.5, 1))
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
             "beta6.p",
             "beta7.p",
             "beta8.p",
             "beta9.p",
             "beta10.p",
             "N",
             "fit",
             "fit.new")

dwri.p.data2 <- list(C = as.matrix(DWRI[, 1:6]), 
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
                      #site = as.numeric(site),
                      Temp.s = Temp.s[ ,1:6],
                      Temp.s2 = Temp.s*Temp.s,
                      RH.s = RH.s[ ,1:6],
                      Precip.s = Precip.s[ ,1:6],
                      VegDamp = VDamp,
                      VegWet = VWet,
                      Mist = Mist,
                      Rain = Rain)

saved.per.chain <- 2000
nc <- 3
nb <- 300000
nt <- 10
ni <- saved.per.chain*nt + nb

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(3) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.p.data2", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_p2.txt", dwri.p.data2, inits, n.adapt = 50000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 1000, thin = 1)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.list.p2 <- mcmc.list(out)
#dwri.list.od2 <- mcmc.list(out[[1]], out[[2]], out[[4]])
plot(dwri.list.p2)
par(mfrow = c(1,1))
summary(dwri.list.p2)

stopCluster(CL)

library(R2jags)
dwri.p2 <- jags(data = dwri.p.data2, parameters.to.save = params, inits = inits, model.file = "dwri_p2.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)

# Posterior Predictive Checks
mean(dwri.p1$BUGSoutput$sims.list$fit.new > dwri.p1$BUGSoutput$sims.list$fit, na.rm = TRUE)
mean(dwri.p1$BUGSoutput$mean$fit) / mean(dwri.p1$BUGSoutput$mean$fit.new)
plot(dwri.p1$BUGSoutput$sims.list$fit, dwri.p1$BUGSoutput$sims.list$fit.new) # model doesn't fit
abline(0, 1, col = 'red')
**plot(dwri.p1$BUGSoutput$sims.list$alpha.lam, dwri.p1$BUGSoutput$sims.list$alpha.p)

traceplot(dwri.p1) # 
par(mfrow = c(1,1))

(dwri.p1.tab <- data.frame(dwri.p1$BUGSoutput$summary[c("alpha.lam",
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
                                                        #"sigma.site",
                                                        "alpha.p",
                                                        "beta1.p",
                                                        "beta2.p",
                                                        "beta3.p",
                                                        "beta4.p",
                                                        "beta5.p",
                                                        # "beta6.p",
                                                        #"beta7.p",
                                                        #"beta8.p",
                                                        #"beta9.p",
                                                        "beta10.p"),
                                                      c(1,2,3,5,7,8,9)]))


#---------------DWRI Reduced overdispersion model-not run-------------
# Define model
sink("dwri_od1.txt")
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
    #beta10.lam ~ dnorm(0, 0.01)
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta10.lam*gcover[i]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -1, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
      # beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0.2, 1),#,
       sigma.p = runif(1, 0.2, 1))#,
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
             #"beta10.lam",
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
             #"delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             D"fit",
             "fit.new")#,
             #"eval",
             #"y.new")

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
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
                  Temp.s = Temp.s[ ,1:6],
                  Temp.s2 = Temp.s*Temp.s,
                  RH.s = RH.s[ ,1:6],
                  Precip.s = Precip.s[ ,1:6])

saved.per.chain <- 2000
nc <- 3
nb <- 100000
nt <- 5
ni <- saved.per.chain*nt + nb

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(4) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)

system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od1.txt", dwri.od.data, inits, n.adapt = 50000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.list <- mcmc.list(out)
plot(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")]) # seems good, maybe slight autocorr
par(mfrow = c(1,1))
summary(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")])

stopCluster(CL)

plot(out[[1]]) # good, reasonable estimate
plot(out[[2]]) # major problems
plot(out[[3]]) # major problem - remove
plot(out[[4]]) # mostly okay but unsure - similar to 1

(out[[1]]) # good, reasonable estimate
plot(out[[2]]) # major problems
plot(out[[3]]) # major problem - remove
plot(out[[4]]) # mostly okay but unsure - similar to 1

CL <- makeCluster(4) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out2 <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od1.txt", dwri.od.data, inits, n.adapt = 100000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 
stopCluster(CL)

plot(out2[[1]]) # good, reasonable estimate
plot(out2[[2]]) # seems reasonable
plot(out2[[3]]) # major problem with N unbounded plus wandering - remove
plot(out2[[4]]) # fine
par(mfrow = c(1,1))

summary(out[[1]])[c("beta1.lam", "beta2.lam")] # no convergence
summary(out[[4]])[c("beta1.lam", "beta2.lam")] # no convergence
summary(out[[2]]) # no convergence

dwri.list <- mcmc.list(out[[1]], out[[4]], out2[[1]], out2[[2]], out2[[4]])
plot(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")]) # seems good, maybe slight autocorr
par(mfrow = c(1,1))
summary(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")])



dwri.od1 <- jags(data = dwri.od.data, parameters.to.save = params, inits = inits, model.file = "dwri_od1.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)

# Posterior Predictive Checks
mean(dwri.od1$BUGSoutput$sims.list$fit.new > dwri.od1$BUGSoutput$sims.list$fit, na.rm = TRUE)
mean(dwri.od1$BUGSoutput$mean$fit) / mean(dwri.od1$BUGSoutput$mean$fit.new)
plot(dwri.od1$BUGSoutput$sims.list$fit, dwri.od1$BUGSoutput$sims.list$fit.new) # excellent fit
abline(0, 1, col = 'red')
plot(dwri.od1$BUGSoutput$sims.list$alpha.lam, dwri.od1$BUGSoutput$sims.list$alpha.p)

plot(dwri.od1$BUGSoutput$sims.list$y.new[1:50, , ], dwri.od1$BUGSoutput$sims.list$eval[1:50, , ])
abline(0, 1, col = 'red')

traceplot(dwri.od1) # 
dev.off()

print(dwri.od1$BUGSoutput$summary[1:50, c(1,2,3, 5, 7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)
# good fit and convergence but large CI making little sig diff from zero (only slope is significant)


# Check for autocorrelation in MCMC saved values
acf(dwri.od1$BUGSoutput$sims.list$alpha.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta1.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta2.lam)
acf(dwri.od1$BUGSoutput$sims.list$alpha.p)
acf(dwri.od1$BUGSoutput$sims.list$beta1.p)
acf(dwri.od1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,50])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,10,6])


sum(dwri.od1$BUGSoutput$mean$N>=1) ## Number Occupied Sites (based on mean)
sum(dwri.od1$BUGSoutput$median$N>=1) ## Number Occupied Sites (based on median)
sum(dwri.od1$BUGSoutput$summary[1:195, 7] >=1) ## Number Occupied Sites (based on 95% CI)
sum(apply(DWRI, 1, sum, na.rm = TRUE) > 0)    # Apparent (naive) distribution among 195 visited sites

# Make table
(SumTab.dwri <- data.frame(dwri.od1$BUGSoutput$summary[c("alpha.lam",
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
                                                         "sigma.site",
                                                      "alpha.p",
                                                      "beta1.p",
                                                      "beta2.p",
                                                      "beta3.p",
                                                      "beta4.p",
                                                      "beta5.p",
                                                      "beta10.p",
                                                         "sigma.p"),
                                                    c(1,2,3,5,7,8,9)]))

log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

Variable <- c("Lambda Intercept", "Elevation", "Elevation2", "Slope", "AspectN", "AspectE", "TPI", "Trail", "LogTWI", "Canopy", "Ground Cover Lambda", "Litter Depth", "Slope2", "Log Stream Dist", "sigma.site", "p Intercept", "Temp", "Temp2", "Precip", "Ground Cover", "Ground Cover2", "RH", "sigma.p")

rownames(SumTab.dwri) <- Variable

print(as.matrix(SumTab.dwri), dig = 2)

write.table(x=SumTab.dwri, file='SummaryTable.dwri.csv', sep=',', row.names=TRUE, col.names=TRUE)


#---------------DWRI full model with overdispersion---not run---------------
# Define model
sink("dwri_od2.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.001)
    beta1.lam ~ dnorm(0, 0.001)
    beta2.lam ~ dnorm(0, 0.001)
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
    
   # for(i in 1:n.sites){
    #eps.lam[i] ~ dnorm(0, tau.site)
    #}
    
    #ls.site ~ dunif(-10, 10)
    #sigma.site <- exp(ls.site)
    #tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.001)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta6.p ~ dnorm(0, 0.01)
    beta7.p ~ dnorm(0, 0.01)
    beta8.p ~ dnorm(0, 0.01)
    beta9.p ~ dnorm(0, 0.001)
    beta10.p ~ dnorm(0, 0.01)
    
    #ls.p ~ dunif(-10, 10)
    #sigma.p <- exp(ls.p)
    #sigma.p ~ dgamma(0.001, 0.001)
    #tau.p <- pow(sigma.p, -2)
    
    #for(i in 1:n.transects){
    #for(j in 1:n.surveys){
    #delta.p[i,j] ~ dnorm(0, tau.p)
    #}
    #}
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta10.lam*gcover[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta6.p*VegDamp[i,j] + beta7.p*VegWet[i,j] + beta8.p*Mist[i,j] + beta9.p*Rain[i,j] + beta10.p*RH.s[i,j] #+ delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
      # beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       #eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0, 1),
       beta6.p = rnorm(1, 0, 1),
       beta7.p = rnorm(1, 0, 1),
       beta8.p = rnorm(1, 0, 1),
       beta9.p = rnorm(1, 0, 1),
       beta10.p = rnorm(1, 0, 1),
       ls.site = rnorm(1, 0, 1))#,#,
       #ls.p = rnorm(1, 0, 1))#,
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
            # "beta10.lam",
             "beta11.lam",
             "beta12.lam",
             "beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta6.p",
             "beta7.p",
             "beta8.p",
             "beta9.p",
             "beta10.p",
             #"eps.lam",
             #"delta.p",
             "sigma.site",
             #"sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data2 <- list(C = as.matrix(DWRI[, 1:6]), 
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6],
                      VegDamp = VDamp,
                      VegWet = VWet,
                      Mist = Mist,
                      Rain = Rain)

saved.per.chain <- 2000
nc <- 3
nb <- 100000
nt <- 5
ni <- saved.per.chain*nt + nb

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(3) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data2", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od2.txt", dwri.od.data2, inits, n.adapt = 10000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 1000, thin = 1)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.list.od2 <- mcmc.list(out)
#dwri.list.od2 <- mcmc.list(out[[1]], out[[2]], out[[4]])
plot(dwri.list.od2)
par(mfrow = c(1,1))
summary(dwri.list.od2)

stopCluster(CL)

plot(out[[1]]) # good, reasonable estimate
plot(out[[2]]) # major problems
plot(out[[3]]) # major problem - remove
plot(out[[4]]) # mostly okay but unsure - similar to 1

(out[[1]]) # good, reasonable estimate
plot(out[[2]]) # major problems
plot(out[[3]]) # major problem - remove
plot(out[[4]]) # mostly okay but unsure - similar to 1

CL <- makeCluster(3) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out2 <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od1.txt", dwri.od.data, inits, n.adapt = 100000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 
stopCluster(CL)

plot(out2[[1]]) # good, reasonable estimate
plot(out2[[2]]) # seems reasonable
plot(out2[[3]]) # major problem with N unbounded plus wandering - remove
plot(out2[[4]]) # fine
par(mfrow = c(1,1))

summary(out[[1]])[c("beta1.lam", "beta2.lam")] # no convergence
summary(out[[4]])[c("beta1.lam", "beta2.lam")] # no convergence
summary(out[[2]]) # no convergence

dwri.list <- mcmc.list(out[[1]], out[[4]], out2[[1]], out2[[2]], out2[[4]])
plot(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # seems good, maybe slight autocorr
par(mfrow = c(1,1))
summary(dwri.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p", "fit", "fit.new")])



dwri.od2 <- jags(data = dwri.od.data, parameters.to.save = params, inits = inits, model.file = "dwri_od1.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)

# Posterior Predictive Checks
mean(dwri.od1$BUGSoutput$sims.list$fit.new > dwri.od1$BUGSoutput$sims.list$fit, na.rm = TRUE)
mean(dwri.od1$BUGSoutput$mean$fit) / mean(dwri.od1$BUGSoutput$mean$fit.new)
plot(dwri.od1$BUGSoutput$sims.list$fit, dwri.od1$BUGSoutput$sims.list$fit.new) # excellent fit
abline(0, 1, col = 'red')
plot(dwri.od1$BUGSoutput$sims.list$alpha.lam, dwri.od1$BUGSoutput$sims.list$alpha.p)

plot(dwri.od1$BUGSoutput$sims.list$y.new[1:50, , ], dwri.od1$BUGSoutput$sims.list$eval[1:50, , ])
abline(0, 1, col = 'red')

traceplot(dwri.od1) # 
dev.off()

print(dwri.od1$BUGSoutput$summary[1:50, c(1,2,3, 5, 7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)
# good fit and convergence but large CI making little sig diff from zero (only slope is significant)


# Check for autocorrelation in MCMC saved values
acf(dwri.od1$BUGSoutput$sims.list$alpha.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta1.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta2.lam)
acf(dwri.od1$BUGSoutput$sims.list$alpha.p)
acf(dwri.od1$BUGSoutput$sims.list$beta1.p)
acf(dwri.od1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,50])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,10,6])


sum(dwri.od1$BUGSoutput$mean$N>=1) ## Number Occupied Sites (based on mean)
sum(dwri.od1$BUGSoutput$median$N>=1) ## Number Occupied Sites (based on median)
sum(dwri.od1$BUGSoutput$summary[1:195, 7] >=1) ## Number Occupied Sites (based on 95% CI)
sum(apply(DWRI, 1, sum, na.rm = TRUE) > 0)    # Apparent (naive) distribution among 195 visited sites

# Make table
(SumTab.dwri <- data.frame(dwri.od1$BUGSoutput$summary[c("alpha.lam",
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
                                                         "sigma.site",
                                                         "alpha.p",
                                                         "beta1.p",
                                                         "beta2.p",
                                                         "beta3.p",
                                                         "beta4.p",
                                                         "beta5.p",
                                                         "beta10.p",
                                                         "sigma.p"),
                                                       c(1,2,3,5,7,8,9)]))

log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

Variable <- c("Lambda Intercept", "Elevation", "Elevation2", "Slope", "AspectN", "AspectE", "TPI", "Trail", "LogTWI", "Canopy", "Ground Cover Lambda", "Litter Depth", "Slope2", "Log Stream Dist", "sigma.site", "p Intercept", "Temp", "Temp2", "Precip", "Ground Cover", "Ground Cover2", "RH", "sigma.p")

rownames(SumTab.dwri) <- Variable

print(as.matrix(SumTab.dwri), dig = 2)

write.table(x=SumTab.dwri, file='SummaryTable.dwri.csv', sep=',', row.names=TRUE, col.names=TRUE)


#---------------DWRI simple model with overdispersion---not run---------------
# Define model
sink("dwri_od3.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.001)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    #beta13.lam ~ dnorm(0, 0.01)
    
     for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    ls.site ~ dunif(-10, 10)
    sigma.site <- exp(ls.site)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta6.p ~ dnorm(0, 0.01)
    beta7.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    ls.p ~ dunif(-10, 10)
    sigma.p <- exp(ls.p)
    #sigma.p ~ dgamma(0.001, 0.001)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta6.p*VegDamp[i,j] + beta7.p*VegWet[i,j] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE))*2 + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -1, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, -0.2, 1),
       beta11.lam = rnorm(1, 0.3, 1),
       #beta13.lam = rnorm(1, 0, 1),
       #eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -2.5, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, -0.2, 1),
       beta3.p = rnorm(1, 0.3, 1),
       beta4.p = rnorm(1, 0.7, 1),
       beta5.p = rnorm(1, -0.2, 1),
       beta6.p = rnorm(1, 1, 1),
       beta7.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 1, 1),
       ls.site = rnorm(1, 0, 1))#,#,
  #ls.p = rnorm(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
             "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             # "beta10.lam",
             "beta11.lam",
             #"beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta5.p",
             "beta6.p",
             "beta7.p",
             "beta10.p",
             #"eps.lam",
             #"delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data3 <- list(C = as.matrix(DWRI[, 1:6]), 
                      n.transects = n.transects, 
                      n.surveys = n.surveys,
                      n.sites = length(site.inits), 
                      elev = elev, 
                      elev2 = elev*elev, 
                      slope = slope,
                      slope2 = slope2, 
                      ltwi = ltwi,
                      trail = trail,
                      canopy = canopy,
                      gcover = gcover,
                      gcover2 = gcover2,
                      litterdepth = litterdepth,
                      #lstream = lstream,
                      site = as.numeric(site),
                      Temp.s = Temp.s[ ,1:6],
                      Temp.s2 = Temp.s*Temp.s,
                      RH.s = RH.s[ ,1:6],
                      Precip.s = Precip.s[ ,1:6],
                      VegDamp = VDamp,
                      VegWet = VWet)

saved.per.chain <- 2000
nc <- 3
nb <- 100000
nt <- 5
ni <- saved.per.chain*nt + nb

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(3) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data3", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od3.txt", dwri.od.data3, inits, n.adapt = 300000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 10000, thin = 10)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.list.od3 <- mcmc.list(out, out2)
plot(dwri.list.od3[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta7.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta6.p", "beta7.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])# looks decent if not great
par(mfrow = c(1,1))
summary(dwri.list.od3[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta7.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta6.p", "beta7.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])#looks good but not convergence for all parameters

stopCluster(CL)

# Posterior Predictive Checks
for(i in 1:3) bayesP <- mean(dwri.list.od3[,"fit.new",][[i]]>dwri.list.od3[,"fit",][[i]]) # Bayesian P value: good fit (0.451)

foo <- as.matrix(dwri.list.od3) # convert to matrix for easier use
plot(foo[, "fit"], foo[,"fit.new"])
abline(0, 1, col = 'red') # excellent fit
par(mfrow = c(1,1))
plot(as.numeric(dwri.list.od3[,"fit.new",]), dwri.list.od3[,"fit",])

traceplot(dwri.list.od3) # 
dev.off()

print(dwri.od1$BUGSoutput$summary[1:50, c(1,2,3, 5, 7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)
# good fit and convergence but large CI making little sig diff from zero (only slope is significant)


# Check for autocorrelation in MCMC saved values
acf(dwri.od1$BUGSoutput$sims.list$alpha.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta1.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta2.lam)
acf(dwri.od1$BUGSoutput$sims.list$alpha.p)
acf(dwri.od1$BUGSoutput$sims.list$beta1.p)
acf(dwri.od1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,50])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,10,6])


sum(dwri.od1$BUGSoutput$mean$N>=1) ## Number Occupied Sites (based on mean)
sum(dwri.od1$BUGSoutput$median$N>=1) ## Number Occupied Sites (based on median)
sum(dwri.od1$BUGSoutput$summary[1:195, 7] >=1) ## Number Occupied Sites (based on 95% CI)
sum(apply(DWRI, 1, sum, na.rm = TRUE) > 0)    # Apparent (naive) distribution among 195 visited sites

# Make table
(SumTab.dwri <- data.frame(dwri.od1$BUGSoutput$summary[c("alpha.lam",
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
                                                         "sigma.site",
                                                         "alpha.p",
                                                         "beta1.p",
                                                         "beta2.p",
                                                         "beta3.p",
                                                         "beta4.p",
                                                         "beta5.p",
                                                         "beta10.p",
                                                         "sigma.p"),
                                                       c(1,2,3,5,7,8,9)]))

log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

Variable <- c("Lambda Intercept", "Elevation", "Elevation2", "Slope", "AspectN", "AspectE", "TPI", "Trail", "LogTWI", "Canopy", "Ground Cover Lambda", "Litter Depth", "Slope2", "Log Stream Dist", "sigma.site", "p Intercept", "Temp", "Temp2", "Precip", "Ground Cover", "Ground Cover2", "RH", "sigma.p")

rownames(SumTab.dwri) <- Variable

print(as.matrix(SumTab.dwri), dig = 2)

write.table(x=SumTab.dwri, file='SummaryTable.dwri.csv', sep=',', row.names=TRUE, col.names=TRUE)



#---------------DWRI simple model with overdispersion--not run----------------
# Define model
sink("dwri_od4.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.001)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    #beta9.lam ~ dnorm(0, 0.01)
    #beta11.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    ls.site ~ dunif(-10, 10)
    sigma.site <- exp(ls.site)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
    #beta6.p ~ dnorm(0, 0.01)
    #beta7.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    #beta11.p ~ dnorm(0, 0.01)
    
    ls.p ~ dunif(-10, 10)
    sigma.p <- exp(ls.p)
    #sigma.p ~ dgamma(0.001, 0.001)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE))*2 + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -2, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       #beta9.lam = rnorm(1, -0.2, 1),
       #beta11.lam = rnorm(1, 0.3, 1),
       #beta13.lam = rnorm(1, 0, 1),
       #eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -2.5, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, -0.2, 1),
       #beta11.p = rnorm(1, 0, 1),
       beta4.p = rnorm(1, 0.7, 1),
      # beta5.p = rnorm(1, -0.2, 1),
       #beta6.p = rnorm(1, 1, 1),
       #beta7.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 1, 1),
       ls.site = rnorm(1, 0, 1))#,#,
  #ls.p = rnorm(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
             "beta7.lam",
             "beta8.lam",
             #"beta9.lam",
             # "beta10.lam",
             #"beta11.lam",
             #"beta13.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta4.p",
             #"beta5.p",
            # "beta6.p",
            # "beta7.p",
             "beta10.p",
            # "beta11.p",
             #"eps.lam",
             #"delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data4 <- list(C = as.matrix(DWRI[, 1:6]), 
                      n.transects = n.transects, 
                      n.surveys = n.surveys,
                      n.sites = length(site.inits), 
                      elev = elev, 
                      elev2 = elev*elev, 
                      ltwi = ltwi,
                      trail = trail,
                      #canopy = canopy,
                      gcover = gcover,
                      #gcover2 = gcover2,
                      #litterdepth = litterdepth,
                      #lstream = lstream,
                      site = as.numeric(site),
                      Temp.s = Temp.s[ ,1:6],
                      Temp.s2 = Temp.s*Temp.s,
                      RH.s = RH.s[ ,1:6])
                      #Precip.s = Precip.s[ ,1:6],
                      #VegDamp = VDamp,
                      #VegWet = VWet)

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(4) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.od.data4", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od4.txt", dwri.od.data4, inits, n.adapt = 300000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 10000, thin = 10)
    return(as.mcmc(fm))
  }) 
}) # 11420s

# Results
dwri.list.od4 <- mcmc.list(out)
plot(dwri.list.od4) # looks good
par(mfrow = c(1,1))
summary(dwri.list.od4) #looks good but not convergence for all parameters

stopCluster(CL)

# Posterior Predictive Checks
for(i in 1:3) bayesP <- mean(dwri.list.od3[,"fit.new",][[i]]>dwri.list.od3[,"fit",][[i]]) # Bayesian P value: good fit (0.451)

foo <- as.matrix(dwri.list.od3) # convert to matrix for easier use
plot(foo[, "fit"], foo[,"fit.new"])
abline(0, 1, col = 'red') # excellent fit
par(mfrow = c(1,1))
plot(as.numeric(dwri.list.od3[,"fit.new",]), dwri.list.od3[,"fit",])

traceplot(dwri.list.od3) # 
dev.off()

print(dwri.od1$BUGSoutput$summary[1:50, c(1,2,3, 5, 7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(dwri.od1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(dwri.od1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)
# good fit and convergence but large CI making little sig diff from zero (only slope is significant)


# Check for autocorrelation in MCMC saved values
acf(dwri.od1$BUGSoutput$sims.list$alpha.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta1.lam)
acf(dwri.od1$BUGSoutput$sims.list$beta2.lam)
acf(dwri.od1$BUGSoutput$sims.list$alpha.p)
acf(dwri.od1$BUGSoutput$sims.list$beta1.p)
acf(dwri.od1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,1])
acf(dwri.od1$BUGSoutput$sims.list$N[ ,50])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(dwri.od1$BUGSoutput$sims.list$delta.p[ ,10,6])


sum(dwri.od1$BUGSoutput$mean$N>=1) ## Number Occupied Sites (based on mean)
sum(dwri.od1$BUGSoutput$median$N>=1) ## Number Occupied Sites (based on median)
sum(dwri.od1$BUGSoutput$summary[1:195, 7] >=1) ## Number Occupied Sites (based on 95% CI)
sum(apply(DWRI, 1, sum, na.rm = TRUE) > 0)    # Apparent (naive) distribution among 195 visited sites

# Make table
(SumTab.dwri <- data.frame(dwri.od1$BUGSoutput$summary[c("alpha.lam",
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
                                                         "sigma.site",
                                                         "alpha.p",
                                                         "beta1.p",
                                                         "beta2.p",
                                                         "beta3.p",
                                                         "beta4.p",
                                                         "beta5.p",
                                                         "beta10.p",
                                                         "sigma.p"),
                                                       c(1,2,3,5,7,8,9)]))

log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

Variable <- c("Lambda Intercept", "Elevation", "Elevation2", "Slope", "AspectN", "AspectE", "TPI", "Trail", "LogTWI", "Canopy", "Ground Cover Lambda", "Litter Depth", "Slope2", "Log Stream Dist", "sigma.site", "p Intercept", "Temp", "Temp2", "Precip", "Ground Cover", "Ground Cover2", "RH", "sigma.p")

rownames(SumTab.dwri) <- Variable

print(as.matrix(SumTab.dwri), dig = 2)

write.table(x=SumTab.dwri, file='SummaryTable.dwri.csv', sep=',', row.names=TRUE, col.names=TRUE)


#---------------DWRI with site effects and no elev2---not run---------------
# Define model
sink("dwri_od5.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.001)
    beta1.lam ~ dnorm(0, 0.001)
    #beta2.lam ~ dnorm(0, 0.001) #can't get convergence with elev2 for DWRI
    beta3.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    
     for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    ls.site ~ dunif(-10, 10)
    sigma.site <- exp(ls.site)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.001)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta6.p ~ dnorm(0, 0.01)
    beta7.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    #ls.p ~ dunif(-10, 10)
    #sigma.p <- exp(ls.p)
    #sigma.p ~ dgamma(0.001, 0.001)
    #tau.p <- pow(sigma.p, -2)
    
    #for(i in 1:n.transects){
    #for(j in 1:n.surveys){
    #delta.p[i,j] ~ dnorm(0, tau.p)
    #}
    #}
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta6.p*VegDamp[i,j] + beta7.p*VegWet[i,j] + beta10.p*RH.s[i,j] #+ delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0, 1),
       beta1.lam = rnorm(1, 2.5, 1),
      # beta2.lam = rnorm(1, -2, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       #eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta6.p = rnorm(1, 0, 1),
       beta7.p = rnorm(1, 0, 1),
       #beta8.p = rnorm(1, 0, 1),
       #beta9.p = rnorm(1, 0, 1),
       beta10.p = rnorm(1, 0, 1),
       ls.site = rnorm(1, 0, 1))#,#,
  #ls.p = rnorm(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             #"beta2.lam", 
             "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             "alpha.p", 
             "beta1.p",
             "beta2.p",
             "beta3.p",
             "beta4.p",
             "beta6.p",
             "beta7.p",
             #"beta8.p",
             #"beta9.p",
             "beta10.p",
             #"eps.lam",
             #"delta.p",
             "sigma.site",
             #"sigma.p",
             "N",
             "p",
             "fit",
             "fit.new")

dwri.data5 <- list(C = as.matrix(DWRI[, 1:6]), 
                      n.transects = n.transects, 
                      n.surveys = n.surveys,
                      n.sites = length(site.inits), 
                      elev = elev, 
                      #elev2 = elev2, 
                      ltwi = ltwi,
                      #tpi = tpi,
                      trail = trail,
                      canopy = canopy,
                      gcover = gcover,
                      #gcover2 = gcover2,
                      litterdepth = litterdepth,
                      lstream = lstream,
                      site = as.numeric(site),
                      Temp.s = Temp.s[ ,1:6],
                      Temp.s2 = Temp.s*Temp.s,
                      RH.s = RH.s[ ,1:6],
                      Precip.s = Precip.s[ ,1:6],
                      VegDamp = VDamp,
                      VegWet = VWet)

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(6) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.data5", "params", "inits", "nc", "ni", "nb", "nt", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_od5.txt", dwri.data5, inits, n.adapt = 1000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 1000, thin = 1)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.list.od5 <- mcmc.list(out)
#dwri.list.od2 <- mcmc.list(out[[1]], out[[2]], out[[4]])
plot(dwri.list.od5[,c("alpha.lam", "beta1.lam", "beta7.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta6.p", "beta7.p", "beta10.p", "sigma.site", "fit", "fit.new")]) 
par(mfrow = c(1,1))
summary(dwri.list.od5[,c("alpha.lam", "beta1.lam", "beta7.lam", "beta8.lam", "beta9.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta6.p", "beta7.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(CL)

#---------------DWRI Reduced w/site effect but no overdispersion model-----not run---------
# Define model
sink("dwri_p6.txt")
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
    #beta10.lam ~ dnorm(0, 0.01)
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
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta10.lam*gcover[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -1, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1))
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
             #"beta10.lam",
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
             #"delta.p",
             "sigma.site",
             #"sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")#,
#"eval",
#"y.new")

dwri.data6 <- list(C = as.matrix(DWRI[, 1:6]), 
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

runif(1)

#out1 <- bugs(data = jags.data, parameters.to.save = params, inits=inits, model.file = "PLJOm1p.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, debug = TRUE)
library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(4) # set nc (usually 3) cores
clusterExport(cl = CL, list("dwri.data6", "params", "inits", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)

system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dwri_p6.txt", dwri.data6, inits, n.adapt = 50000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dwri.p6.list <- mcmc.list(out)
plot(dwri.p6.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site")]) #
par(mfrow = c(1,1))
summary(dwri.p6.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site")])

stopCluster(CL)

#----------DWRI - Overdispersion in abundance and detection linear elev - good model but need to rerun with sigma = 10 not 3 for a prior---not great fit----------
# Define model
sink("dwri_od3.txt")
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
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 10)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
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
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta2.lam*elev2[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j] <- p[i,j] * N[i]
    Residuals[i,j] <- C[i,j] - eval[i,j]
    E[i,j] <- pow((Residuals[i,j]),2) / (eval[i,j] + 0.5)
    
    # Generate replicate data and compute fit stats for them
    y.new[i,j] ~ dbin(p[i,j], N[i])
    Residuals.new[i,j] <- y.new[i,j] - eval[i,j]
    E.new[i,j] <- pow((Residuals.new[i,j]),2) / (eval[i,j] + 0.5) 
    
    diff.cy[i,j] <- C[i,j] - y.new[i,j]
    
    # Squared residuals for SSE
    Res2[i,j] <- pow(Residuals[i,j], 2)
    Res2.new[i,j] <- pow(Residuals.new[i,j], 2)

    # Freeman-Tukey Fit Test
    FT[i,j] <- pow(C[i,j],2) - pow(eval[i,j], 2)
    FT.new[i,j] <- pow(y.new[i,j], 2) - pow(eval[i,j], 2)
    }
    }
    
    # Derived quantities
    totalN<- sum(N[])
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])
    
    SSE <- sum(Res2[,])
    SSE.new <- sum(Res2.new[,])

    FreeTukey <- sum(FT[,])
    FreeTukey.new <- sum(FT.new[,])
    }
    ", fill = TRUE)
sink()

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

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
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
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_od3.txt", dwri.od.data, inits, n.adapt=3000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=10000, thin=4) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_od3 <- mcmc.list(out)
plot(dwri_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.dwri3 <- mean(dwri_od3[, "fit.new",][[i]] > dwri_od3[, "fit",][[i]]) # 0.456 excellent
print(bayesP.dwri3, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(dwri_od3[, "fit",]), as.matrix(dwri_od3[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=dwri_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # not converged even though chains look good

#----------DWRI - simple model w/ linear elev ---good model- USE This Model--------------
# Define model
sink("dwri_od4.txt")
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
    
    sigma.site ~ dunif(0, 10)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 10)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta2.lam*elev2[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
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
       #beta5.p = rnorm(1, 0.5, 1),
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
             #"beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
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
                     #gcover2 = gcover2,
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
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_od4.txt", dwri.od.data, inits, n.adapt=5000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=2000, thin=1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_od4 <- mcmc.list(out)
plot(dwri_od4[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od4[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

save(dwri_od4, file = "Results/JAGS/dwri_mcmc_out.RData")
saveRDS(dwri_od4, file = "Results/JAGS/dwri_mcmc_out.rds")
save.image("data_and_results.RData")

# Check fit
for(i in 1:3) bayesP.dwri4 <- mean(dwri_od4[, "fit.new",][[i]] > dwri_od4[, "fit",][[i]]) # 0.454 excellent
print(bayesP.dwri4, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(dwri_od4[, "fit",]), as.matrix(dwri_od4[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=dwri_od4[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # good convergence

Quants.dwri <- apply(as.matrix(dwri_od4[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.dwri <- apply(as.matrix(dwri_od4[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.dwri <- apply(as.matrix(dwri_od4[ , c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

Dwri.variables <- c("N-intercept", "Elevation", "Slope", "TPI", "Trail", "Log(TWI)", "Canopy", "Ground Cover", "Litter Depth", "Log(Stream Dist)", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Rel. Humidity", "Site SD", "Fit-Data", "Fit-Ideal")

Dwri.summary <- data.frame(Dwri.variables, Means.dwri, SDs.dwri, Quants.dwri["2.5%", ], Quants.dwri["50%", ], Quants.dwri["97.5%", ])

colnames(Dwri.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(Dwri.summary, file = "Dwri_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.dwri <- matrix(NA, 195, 3)
for(i in 1:195){
  N.dwri[i, ] <- apply(as.matrix(dwri_od4[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
}

colnames(N.dwri) <- c("CI_2.5", "Median", "CI_97.5")
N.dwri
cbind(NaiveRichness, DWRImin, N.Dwri)


#----------DWRI - simple model w/ linear elev & hierarchical centering -----------------
# Define model
sink("dwri_od5.txt")
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
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 10)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
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
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta2.lam*elev2[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
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
       eps.lam = rnorm(length(unique(site)), 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       #beta5.p = rnorm(1, 0.5, 1),
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
             #"beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
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
                     #gcover2 = gcover2,
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
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_od5.txt", dwri.od.data, inits, n.adapt=1000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=5000, thin=4) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_od5 <- mcmc.list(out)
plot(dwri_od5[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od5[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)



#----------DWRI - simple model w/ linear elev & hierarchical centering -----------------
# Define model
sink("dwri_od6.txt")
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
    eps.lam[i] ~ dnorm(alpha.lam, tau.site)
    }
    
    sigma.site ~ dunif(0, 10)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
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
    
    log(lambda[i]) <- beta1.lam*elev[i] + beta3.lam*slope[i] + beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta2.lam*elev2[i]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
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

Nst <- apply(DWRI, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 4, 1),
       # beta2.lam = rnorm(1, -2, 1),
       beta3.lam = rnorm(1, 0, 1),
       #beta4.lam = rnorm(1, 0, 1),
       #beta5.lam = rnorm(1, 0, 1),
       # beta6.lam = rnorm(1, 0, 1),
       # beta7.lam = rnorm(1, 0, 1),
       # beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       beta10.lam = rnorm(1, 0, 1),
       # beta11.lam = rnorm(1, 0, 1),
       #beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       #beta5.p = rnorm(1, 0.5, 1),
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
             #"beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")

dwri.od.data <- list(C = as.matrix(DWRI[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(unique(site)), 
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
                     #gcover2 = gcover2,
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
clusterExport(cl, c("dwri.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_od6.txt", dwri.od.data, inits, n.adapt=5000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=2000, thin=1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_od6 <- mcmc.list(out)
plot(dwri_od6[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(dwri_od6[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)




#----------DWRI Occupancy Model----------
# Define model
sink("dwri_od3_occ.txt")
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
    
    sigma.site ~ dunif(0, 25)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
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
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    lpsi.lim[i] <- min(999, max(-999, lpsi[i])) # Help stabilize the logit
    psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))    
    
    for(j in 1:n.surveys){
    y[i, j] ~ dbern(eff.p[i, j])
    eff.p[i,j] <- z[i] * p[i,j]
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j] + delta.p[i,j]
    
    } 
    }
    }
    ", fill = TRUE)
sink()

# Convert DWRI Data to 0/1
DWRI01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    DWRI01[i,j] <- ifelse(DWRI[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(DWRI01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       #beta2.lam = rnorm(1, 0, 1),
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
       #beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1),
       sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
            # "beta2.lam", 
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
             #"beta5.p",
             "beta10.p",
             #"eps.lam",
             "delta.p",
             "sigma.site",
             "sigma.p",
             "psi",
             "z")

dwri.od.data01 <- list(y = as.matrix(DWRI01[, 1:6]), 
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
                      # gcover2 = gcover2,
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
clusterExport(cl, c("dwri.od.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_od3_occ.txt", dwri.od.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_occ_od <- mcmc.list(out)
plot(dwri_occ_od [,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")]) # 
par(mfrow = c(1,1))
summary(dwri_occ_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")])

stopCluster(cl)

print(gelman.diag(x=dwri_occ_od[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")]), dig=3) # 


#----------DWRI Occupancy Model----------
# Define model
sink("dwri_occ.txt")
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
    
    sigma.site ~ dunif(0, 25)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    #beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    # Likelihood
    for(i in 1:n.transects){
    z[i] ~ dbern(psi[i])
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    lpsi.lim[i] <- min(999, max(-999, lpsi[i])) # Help stabilize the logit
    psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))    
    
    for(j in 1:n.surveys){
    y[i, j] ~ dbern(eff.p[i, j])
    eff.p[i,j] <- z[i] * p[i,j]
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta10.p*RH.s[i,j]
    
    } 
    }
    }
    ", fill = TRUE)
sink()

# Convert DWRI Data to 0/1
DWRI01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    DWRI01[i,j] <- ifelse(DWRI[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(DWRI01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       #beta2.lam = rnorm(1, 0, 1),
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
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, 0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       #beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             # "beta2.lam", 
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
             #"beta5.p",
             "beta10.p",
             #"eps.lam",
             "sigma.site",
             "psi",
             "z")

dwri.data01 <- list(y = as.matrix(DWRI01[, 1:6]), 
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
                       # gcover2 = gcover2,
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
clusterExport(cl, c("dwri.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("dwri_occ.txt", dwri.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
dwri_occ <- mcmc.list(out)
plot(dwri_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")]) # 
par(mfrow = c(1,1))
summary(dwri_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")])

stopCluster(cl)

print(gelman.diag(x=dwri_occ[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta10.p", "sigma.site")]), dig=3) # 


#######################################
# Eurycea wilderae - Bayesian Analysis
#######################################

#---------------EWIL full overdispersion model -No convergence -----------------
# Define model
sink("ewil_od1.txt")
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
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

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       beta2.lam = rnorm(1, 0, 1),
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("ewil_od1.txt", ewil.od.data, inits, n.adapt=300000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
ewil_od1 <- mcmc.list(out)
plot(ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.ewil1 <- mean(ewil_od1[, "fit.new",][[i]] > ewil_od1[, "fit",][[i]]) # 
print(bayesP.ewil1, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(ewil_od1[, "fit",]), as.matrix(ewil_od1[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # 

Quants.ewil <- apply(as.matrix(ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.ewil <- apply(as.matrix(ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), 2, FUN = mean)

SDs.ewil <- apply(as.matrix(ewil_od1[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta12.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), 2, FUN = sd)

  
  Ewil.variables <- c("N-intercept", "Elevation", "Elevation^2", "Slope", "Aspect-N", "Aspect-E", "TPI", "Trail", "Log(TWI)", "Canopy", "Ground Cover", "Litter Depth", "Slope^2", "Log(Stream Dist)", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Site SD", "Fit-Data", "Fit-Ideal")
  
  Ewil.summary <- data.frame(Ewil.variables, Means.ewil, SDs.ewil, Quants.ewil["2.5%", ], Quants.ewil["50%", ], Quants.ewil["97.5%", ])
  
  colnames(Ewil.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")
  
  write.table(Ewil.summary, file = "Ewil_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)
  
  N.ewil <- matrix(NA, 195, 3)
  for(i in 1:195){
    N.ewil[i, ] <- apply(as.matrix(ewil_od1[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
  }
  
  colnames(N.ewil) <- c("CI_2.5", "Median", "CI_97.5")
  N.ewil
  cbind(NaiveRichness, EWILmin, N.ewil)
  
  
#---------R2Jags------------
saved.per.chain <- 1000
nc <- 3
nb <- 300000
nt <- 5
ni <- saved.per.chain*nt + nb

runif(1)

ewil.od1 <- jags(data = ewil.od.data, parameters.to.save = params, inits = inits, model.file = "ewil_od1.txt", n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)

# Posterior Predictive Checks
mean(ewil.od1$BUGSoutput$sims.list$fit.new > ewil.od1$BUGSoutput$sims.list$fit, na.rm = TRUE)
mean(ewil.od1$BUGSoutput$mean$fit) / mean(ewil.od1$BUGSoutput$mean$fit.new)
plot(ewil.od1$BUGSoutput$sims.list$fit, ewil.od1$BUGSoutput$sims.list$fit.new) # 
abline(0, 1, col = 'red')
plot(ewil.od1$BUGSoutput$sims.list$alpha.lam, ewil.od1$BUGSoutput$sims.list$alpha.p)

plot(ewil.od1$BUGSoutput$sims.list$y.new[1:50, , ], ewil.od1$BUGSoutput$sims.list$eval[1:50, , ])
abline(0, 1, col = 'red')

traceplot(ewil.od1) # good mixing and convergence
dev.off()

print(ewil.od1$BUGSoutput$summary[1:50, c(1,2,3, 5, 7,8,9)], dig = 3)
print(ewil.od1$BUGSoutput$summary[51:100, c(1,2,3,7,8,9)], dig = 3)
print(ewil.od1$BUGSoutput$summary[1:195, c(1,2,3,7,8,9)], dig = 3)
print(ewil.od1$BUGSoutput$summary[196:393, c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["alpha.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta1.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta2.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta3.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta4.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta5.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta6.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta7.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta8.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta9.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta10.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta11.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta12.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["beta13.lam", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["alpha.p", c(1,2,3,7,8,9)], dig = 2)
print(ewil.od1$BUGSoutput$summary["eval[1,1]", c(1,2,3,7,8,9)], dig = 2)
# good fit and convergence but large CI making little sig diff from zero (only slope is significant)


# Check for autocorrelation in MCMC saved values
acf(ewil.od1$BUGSoutput$sims.list$alpha.lam)
acf(ewil.od1$BUGSoutput$sims.list$beta1.lam)
acf(ewil.od1$BUGSoutput$sims.list$beta2.lam)
acf(ewil.od1$BUGSoutput$sims.list$alpha.p)
acf(ewil.od1$BUGSoutput$sims.list$beta1.p)
acf(ewil.od1$BUGSoutput$sims.list$eps.lam[ ,1])
acf(ewil.od1$BUGSoutput$sims.list$N[ ,1])
acf(ewil.od1$BUGSoutput$sims.list$N[ ,50])
acf(ewil.od1$BUGSoutput$sims.list$delta.p[ ,1,1])
acf(ewil.od1$BUGSoutput$sims.list$delta.p[ ,2,2])
acf(ewil.od1$BUGSoutput$sims.list$delta.p[ ,10,6])


sum(ewil.od1$BUGSoutput$mean$N>=1) ## Number Occupied Sites (based on mean)
sum(ewil.od1$BUGSoutput$median$N>=1) ## Number Occupied Sites (based on median)
sum(ewil.od1$BUGSoutput$summary[1:195, 7] >=1) ## Number Occupied Sites (based on 95% CI)
sum(ewil.od1$BUGSoutput$summary[1:195, 3] >=1) ## Number Occupied Sites (based on lower CI)
sum(apply(EWIL, 1, sum, na.rm = TRUE) > 0)    # Apparent (naive) distribution among 195 visited sites

# Make table
(SumTab.ewil <- data.frame(ewil.od1$BUGSoutput$summary[c("alpha.lam",
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
                                                         "sigma.site",
                                                         "alpha.p",
                                                         "beta1.p",
                                                         "beta2.p",
                                                         "beta3.p",
                                                         "beta4.p",
                                                         "beta5.p",
                                                         "beta10.p",
                                                         "sigma.p"),
                                                       c(1,2,3,5,7,8,9)]))

log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

rownames(SumTab.ewil) <- Variable

print(as.matrix(SumTab.ewil), dig = 2)

write.table(x=SumTab.ewil, file='SummaryTable.ewil.csv', sep=',', row.names=TRUE, col.names=TRUE)


range(Data$elev[which(EWILmin > 0)], na.rm=TRUE) # observed Ewil 447 - 2020 m (Naive)
range(Data$elev[which(ewil.od1$BUGSoutput$median$N > 0)], na.rm=TRUE) # median at sites 447 - 2020 


#-----------Figures-----------------------------------

meanc <- fixef(lmm1) # glmm analysis of count data
lmc <- coef(lm1) # glm analysis


# Eleveation effect on Abundance at mean other conditions measured
Elev <- seq(from = 400, to = max(Data$elev), length.out = 1000)
Elev.s <- (Elev - mean(Data$elev))/sd(Data$elev)
a.lam.ewil <- ewil.od1$BUGSoutput$mean$alpha.lam
b.elev.ewil <- ewil.od1$BUGSoutput$mean$beta1.lam
b.elev2.ewil <- ewil.od1$BUGSoutput$mean$beta2.lam
N.elev.ewil <- exp(a.lam.ewil + b.elev.ewil*Elev.s + b.elev2.ewil*Elev.s^2)
c.elev <- exp(meanc[1] + meanc["elevl"]*Elev.s + meanc["elev2l"]*Elev.s^2)
lm.elev <- exp(lmc[1] + lmc["elev"]*Elev.s + lmc["elev2"]*Elev.s^2)

# Litter effect on Abundance at mean other conditions measured
Lstream <- seq(from = min(Data$lstream, na.rm = TRUE), to = max(Data$lstream, na.rm = TRUE), length.out = 100)
Lstream.s <- (Lstream - mean(Data$lstream))/sd(Data$lstream)
N.lstream <- exp(a.lam.ewil + ewil.od1$BUGSoutput$mean$beta13.lam*Lstream.s)
c.litter <- exp(meanc[1] + meanc["litterl"]*Litter.s)
lm.litter <- exp(lmc[1] + lmc["litterdepth"]*Litter.s)

# PLOTS
par(mfrow = c(1,2), oma=c(2,2,0,0), mar = c(4,3,2,2) + 0.1)
plot(Elev, N.elev.ewil, type = "l",
     ylim = c(0, 7),
     xlab = 'Elevation (m)',
     ylab = '')
#points(Data$elev, PLJOmin)
#lines(Elev, c.elev, col = 'blue', lty = 2) # GLMM prediction without detection
#lines(Elev, lm.elev, col = 'red', lty = 3) # GLM prediction for 1st visit

plot(exp(Lstream), N.lstream, type = "l",
     ylim = c(0, 7),
     xlab = 'Stream distance (m)',
     ylab = '')
#lines(Slo, c.slope, col = 'blue', lty = 2)
#lines(Slo, lm.slope, col = 'red', lty = 3)

mtext(expression(paste("Abundance of ", italic("E. wilderae"))), outer = TRUE, side=2 )
par(mfrow = c(1,1))

# elevation of peak abundance 
Elev[which(N.elev == max(N.elev))] # 1547 m
Elev[which(c.elev == max(c.elev))] # 1530 m
Elev[which(lm.elev == max(lm.elev))] # 1448 m

# min elevation of occurence
max(Elev[which(N.elev < 1)]) # predicted min elevation under average conditions = 826 m

Medians <- out.m2.od2$BUGSoutput$median
range(Data$elev[which(Medians$N == 0)]) # predicted zeros occur from 412 - 946 m (depends on other habitat conditions)

dev.off()

# Don't include canopy plot since not significant (line will just look flat and doesn't technically differ from 0 - not worth taking up space in manuscript - just show beta in table) - although sig for the naive models
plot(Canopy, N.canopy, type = "l",
     ylim = c(0, 75),
     xlab = 'Canopy',
     ylab = expression(paste("Abundance of ", italic("P. jordani"))))
lines(Canopy, c.canopy, col = 'blue', lty = 2)
lines(Canopy, lm.canopy, col = 'red', lty = 3)


# Detection Plots Temperature, precipitation, ground cover, and RH 
library(boot)
Temp <- c(Data$temp1, Data$temp2, Data$temp3, Data$temp4, Data$temp5, Data$temp6)
TempR <- seq(from = min(Temp, na.rm = TRUE), to = max(Temp, na.rm = TRUE), length.out = 100)
Temp.s <- (TempR - mean(unlist(Temp), na.rm=TRUE))/sd(unlist(Temp), na.rm=TRUE)
p.temp <- inv.logit(ewil.od1$BUGSoutput$mean$alpha.p + ewil.od1$BUGSoutput$mean$beta1.p*Temp.s + ewil.od1$BUGSoutput$mean$beta2.p*Temp.s^2)

Precip <- c(Data$precip1, Data$precip2, Data$precip3, Data$precip4, Data$precip5, Data$precip6)
Precipr <- seq(from = min(Precip, na.rm = TRUE), to = max(Precip, na.rm = TRUE), length.out = 100)
Precip.s <- (Precipr - mean(unlist(Precip), na.rm=TRUE))/sd(unlist(Precip), na.rm = TRUE)
p.precip <- inv.logit(ewil.od1$BUGSoutput$mean$alpha.p + ewil.od1$BUGSoutput$mean$beta3.p*Precip.s)

RH <- c(Data$RH1, Data$RH2, Data$RH3, Data$RH4, Data$RH5, Data$RH6)
RHr <- seq(from = min(RH, na.rm = TRUE), to = max(RH, na.rm = TRUE), length.out = 100)
RH.s <- (RHr - mean(unlist(RH), na.rm=TRUE))/sd(unlist(RH), na.rm = TRUE)
p.RH <- inv.logit(ewil.od1$BUGSoutput$mean$alpha.p + ewil.od1$BUGSoutput$mean$beta10.p*RH.s)

p.cover <- inv.logit(ewil.od1$BUGSoutput$mean$alpha.p + ewil.od1$BUGSoutput$mean$beta4.p*Cover.s + ewil.od1$BUGSoutput$mean$beta5.p*Cover.s^2)

par(mfrow = c(2,2), oma=c(2,2,0,0), mar = c(4,3,2,2) + 0.1)
plot(TempR, p.temp, type = "l",
     ylim = c(0, 0.5),
     xlab = 'Temperature (C)',
     ylab = '')

plot(Precipr, p.precip, type = "l",
     ylim = c(0, 0.5),
     xlab = 'Precipitation past 24 hrs (cm)',
     ylab = '')

plot(RHr, p.RH, type = "l",
     ylim = c(0, 0.5),
     xlab = 'Relative Humidity (%)',
     ylab = '')

plot(Cover, p.cover, type = "l",
     ylim = c(0, 0.5),
     xlab = 'Vegetative ground cover (%)',
     ylab = '')

mtext(expression(paste("Detection probability")), outer = TRUE, side=2 )
par(mfrow=c(1,1))

p.temp[which(p.temp == max(p.temp))]
TempR[which(p.temp == max(p.temp))] # Ideal Temperature = 21.2 C results in 31.2% detection probability with all else at mean conditions



#---------------EWIL Reduced model with random site effect on N but no overdispersion------not run-----------
# Define model
sink("ewil_p2.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
    beta4.lam ~ dnorm(0, 0.01)
    beta5.lam ~ dnorm(0, 0.01)
    beta6.lam ~ dnorm(0, 0.01)
    beta7.lam ~ dnorm(0, 0.01)
    beta8.lam ~ dnorm(0, 0.01)
    beta9.lam ~ dnorm(0, 0.01)
    beta11.lam ~ dnorm(0, 0.01)
    beta13.lam ~ dnorm(0, 0.01)
    
    for(i in 1:n.sites){
    eps.lam[i] ~ dnorm(0, tau.site)
    }
    
    ls.site ~ dunif(-10, 10)
    sigma.site <- exp(ls.site)
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
    for(j in 1:n.surveys){
    C[i, j] ~ dbin(p[i, j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))    
    lp.lim[i,j] <- min(999, max(-999, lp[i,j])) # Help stabilize the logit
    
    lp[i, j] <- alpha.p + beta1.p*Temp.s[i,j] + beta2.p*Temp.s2[i,j] + beta3.p*Precip.s[i,j] + beta4.p*gcover[i] + beta5.p*gcover2[i] + beta10.p*RH.s[i,j]
    
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

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 1, 1),
       beta1.lam = rnorm(1, 3, 1),
       beta2.lam = rnorm(1, -3, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
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
       beta10.p = rnorm(1, 0.5, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
            
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "beta7.lam",
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
             "sigma.site",
             "N",
             "p",
             "fit",
             "fit.new")

ewil.data2 <- list(C = as.matrix(EWIL[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
                     #slope = slope,
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])


runif(1)


library(parallel)
# Run JAGS in parallel so each chain in run simulaneously using a separate computer core (CPU)
CL <- makeCluster(7) # set nc (usually 3) cores
clusterExport(cl = CL, list("ewil.data2", "params", "inits", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)
system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("ewil_p2.txt", ewil.data2, inits, n.adapt = 100000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 10000, thin = 10)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
ewil.list2 <- mcmc.list(out)
plot(ewil.list2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # a few bad chains & problems with elev & elev2
par(mfrow = c(1,1))
summary(ewil.list2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) 

stopCluster(CL)

plot(ewil.list2[,c("alpha.lam", "beta1.lam", "beta2.lam")])
plot(out[[1]][,c("alpha.lam", "beta1.lam")]) # get rid of chain 1
plot(out[[2]][,c("alpha.lam", "beta1.lam")]) # bad
plot(out[[3]][,c("alpha.lam", "beta1.lam")])
plot(out[[4]][,c("alpha.lam", "beta1.lam")])
plot(out[[5]][,c("alpha.lam", "beta1.lam")]) # bad
plot(out[[6]][,c("alpha.lam", "beta1.lam")]) # badish
plot(out[[7]][,c("alpha.lam", "beta1.lam")]) # okay

# I guess it's more than just a few problems. Just use overdispersed detection model for inference even though it has high uncertainty.

    
    

    )

#---------------EWIL reduced overdispersion model - Good model-----------------
# Define model
sink("ewil_od2.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
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

Nst <- apply(EWIL, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       beta2.lam = rnorm(1, 0, 1),
       beta3.lam = rnorm(1, 0, 1),
      # beta4.lam = rnorm(1, 0, 1),
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
             "beta2.lam", 
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
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])
library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("ewil_od2.txt", ewil.od.data, inits, n.adapt=300000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
ewil_od2 <- mcmc.list(out)
plot(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # Only "problem" is perfect inverse correlation between elev and elev^2
par(mfrow = c(1,1))
summary(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

save(ewil_od2, file = "Results/JAGS/ewil_mcmc_out.RData")

# Check fit
for(i in 1:3) bayesP.ewil2 <- mean(ewil_od2[, "fit.new",][[i]] > ewil_od2[, "fit",][[i]]) # 0.499 perfect
print(bayesP.ewil2, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(ewil_od2[, "fit",]), as.matrix(ewil_od2[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # All < 1.1 w/95% CI

# This will be very slow - don't rerun
gelman.diag(ewil_od2)

acf(as.matrix(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]))

# Check for excessive autocorrelation in chains
par(mfrow = c(1,1))
acf(as.matrix(ewil_od2[,"alpha.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta1.lam"]), lag.max = 1000) # high autocorr
acf(as.matrix(ewil_od2[,"beta2.lam"]), lag.max = 1000) # High autocorr
acf(as.matrix(ewil_od2[,"beta3.lam"]), lag.max = 1000) # slight autocorr
acf(as.matrix(ewil_od2[,"beta6.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta7.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta8.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta9.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta10.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta11.lam"]), lag.max = 1000)
acf(as.matrix(ewil_od2[,"beta13.lam"]), lag.max = 1000)

# Effective Size
effectiveSize(ewil_od2[,"alpha.lam"])
effectiveSize(ewil_od2[,"beta1.lam"])
effectiveSize(ewil_od2[,"beta2.lam"])
effectiveSize(ewil_od2[,"beta3.lam"])
effectiveSize(ewil_od2[,"beta6.lam"])
effectiveSize(ewil_od2[,"beta7.lam"])
effectiveSize(ewil_od2[,"beta8.lam"])
effectiveSize(ewil_od2[,"beta9.lam"])

# Rejection Rate
rejectionRate(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

# Acceptance Rate
print(1 - rejectionRate(ewil_od2[,c("alpha.lam", "beta1.lam")]), dig = 6)

# Summaries
Quants.ewil <- apply(as.matrix(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.ewil <- apply(as.matrix(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), 2, FUN = mean)

SDs.ewil <- apply(as.matrix(ewil_od2[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta10.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]), 2, FUN = sd)

Ewil.variables <- c("N-intercept", "Elevation", "Elevation^2", "Slope", "TPI", "Trail", "Log(TWI)", "Canopy", "Ground Cover", "Litter Depth", "Log(Stream Dist)", "p-intercept", "Temperature", "Temperature^2", "24-hr Precip", "Ground Cover", "Ground Cover^2", "Rel. Humidity", "Site SD", "Detection SD", "Fit-Data", "Fit-Ideal")
  
Ewil.summary <- data.frame(Ewil.variables, Means.ewil, SDs.ewil, Quants.ewil["2.5%", ], Quants.ewil["50%", ], Quants.ewil["97.5%", ])
  
colnames(Ewil.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")
  
write.table(Ewil.summary, file = "Ewil_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)
  
N.ewil <- matrix(NA, 195, 3)
for(i in 1:195){
    N.ewil[i, ] <- apply(as.matrix(ewil_od2[, c(paste("N[", i, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.025, 0.5, 0.975))
  }
  
  colnames(N.ewil) <- c("CI_2.5", "Median", "CI_97.5")
  N.ewil
  cbind(NaiveRichness, EWILmin, N.ewil)

write.table(cbind(NaiveRichness, EWILmin, N.ewil), file = "Ewil_N.csv", sep = ",", row.names = TRUE)
  
# Peak Elevation
el <- seq(400, 2025, by = 1)
el2 <- el^2
el.s <- (el - mean(Data$elev)) / sd(Data$elev)
el2.s <- (el2 - mean(Data$elev2)) / sd(Data$elev2)
df.el <- data.frame(el)
  
# convert mcmc list into dataframe
ewil.iter <- as.data.frame(as.matrix(ewil_od2))

# filter to elevation columns
ewil.N.elev <- ewil.iter %>%
  dplyr::select(Nintercept = alpha.lam, Elevation = beta1.lam, Elevation2 = beta2.lam) 

# calculate abundance at each iteration at 
ewil.peak.N <- rep(NA_real_, times = nrow(ewil.iter))
for(i in 1:nrow(ewil.iter)) {
  df.el[ , paste0("N",i)] <- exp(ewil.N.elev[i, "Nintercept"] + ewil.N.elev[i, "Elevation"] * el.s + ewil.N.elev[i, "Elevation2"] * el2.s)
  ewil.peak.N[i] <- df.el[which(df.el[ , paste0("N",i)] == max(df.el[ , paste0("N",i)])), "el"]
}

# calculate expected peak elevation with 95% CRI
median(ewil.peak.N)
quantile(ewil.peak.N, probs = c(0.025, 0.975))

# summarize iterations for plotting
df.el <- df.el %>%
  dplyr::select(-el)
df.el.summary <- data.frame(el)
df.el.summary$median <- apply(df.el, MARGIN = 1, FUN = median)
df.el.summary$LCRI <- apply(df.el, MARGIN = 1, FUN = quantile, probs = c(0.025))
df.el.summary$UCRI <- apply(df.el, MARGIN = 1, FUN = quantile, probs = c(0.975))

# plot elevation vs. expected abundance with CRI
ewil.elev.N.plot <- ggplot(df.el.summary, aes(el, median)) + geom_line() + xlab("Elevation (m)") + ylab("Abundance") + ggtitle(expression(italic("E. wilderae"))) + geom_ribbon(aes(ymin=LCRI, ymax=UCRI),alpha=0.3) + theme_bw()


#----------EWIL Occupancy Model--can't handle site effects--------
# Define model
sink("ewil_od1_occ.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    
    sigma.site ~ dunif(0, 10)
    tau.site <- 1/(sigma.site*sigma.site)
    
    alpha.p ~ dnorm(0, 0.01)
    beta1.p ~ dnorm(0, 0.01)
    beta2.p ~ dnorm(0, 0.01)
    beta3.p ~ dnorm(0, 0.01)
    beta4.p ~ dnorm(0, 0.01)
    beta5.p ~ dnorm(0, 0.01)
    beta10.p ~ dnorm(0, 0.01)
    
    sigma.p ~ dunif(0, 10)
    tau.p <- pow(sigma.p, -2)
    
    for(i in 1:n.transects){
    for(j in 1:n.surveys){
    delta.p[i,j] ~ dnorm(0, tau.p)
    }
    }
    
    # Likelihood
    for(i in 1:n.transects){
    z[i] ~ dbern(psi[i])
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]

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
EWIL01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    EWIL01[i,j] <- ifelse(EWIL[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(EWIL01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       beta2.lam = rnorm(1, 0, 1),
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
             "beta2.lam", 
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
             "psi",
             "z",
             #"p",
             "fit",
             "fit.new")

ewil.od.data01 <- list(y = as.matrix(EWIL01[, 1:6]), 
                     n.transects = n.transects, 
                     n.surveys = n.surveys,
                     n.sites = length(site.inits), 
                     elev = elev, 
                     elev2 = elev2, 
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.od.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("ewil_od1_occ.txt", ewil.od.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
ewil_occ_od <- mcmc.list(out)
plot(ewil_occ_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_occ_od[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.ewil_occ_od <- mean(ewil_occ_od[, "fit.new",][[i]] > ewil_occ_od[, "fit",][[i]]) # 
print(bayesP.ewil_occ_od, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(ewil_occ_od[, "fit",]), as.matrix(ewil_occ_od[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=ewil_occ_od[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # 

#----------EWIL Occupancy Model---okay - can't handle random site effects although maybe useable-------
# run without overdispersed detection because diff detection than abund models. It could help reduce uncertainty in detection estimates.
# Define model
sink("ewil_occ.txt")
cat("
    model{
    # Priors
    alpha.lam ~ dnorm(0, 0.01)
    beta1.lam ~ dnorm(0, 0.01)
    beta2.lam ~ dnorm(0, 0.01)
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
    
    sigma.site ~ dunif(0, 10)
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
    
    lpsi[i] <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
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
    fit.new <- sum(Presi.new[,]) # Discrepancy for replicate data set
    }
    ", fill = TRUE)
sink()

# Convert EWIL Data to 0/1
EWIL01 <- matrix(NA, 195, 6)
for(i in 1:195){
  for(j in 1:6){
    EWIL01[i,j] <- ifelse(EWIL[i,j] == 0, 0, 1)
  }
} 

Zst <- apply(EWIL01, 1, function(x) max(x, na.rm = TRUE))
inits <- function(){
  list(z = Zst,
       alpha.lam = rnorm(1, 0.5, 1),
       beta1.lam = rnorm(1, 0, 1),
       beta2.lam = rnorm(1, 0, 1),
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
       sigma.site = runif(1, 0, 1))
      # sigma.p = runif(1, 0, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
             "beta2.lam", 
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
             "sigma.site",
             "psi",
             "z",
             #"p",
             "fit",
             "fit.new")

ewil.data01 <- list(y = as.matrix(EWIL01[, 1:6]), 
                       n.transects = n.transects, 
                       n.surveys = n.surveys,
                       n.sites = length(site.inits), 
                       elev = elev, 
                       elev2 = elev2, 
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
                       Temp.s = Temp.s[ ,1:6],
                       Temp.s2 = Temp.s*Temp.s,
                       RH.s = RH.s[ ,1:6],
                       Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("ewil.data01", "inits", "params", "Zst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("ewil_occ.txt", ewil.data01, inits, n.adapt=50000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
ewil_occ <- mcmc.list(out)
plot(ewil_occ[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(ewil_occ[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")])

stopCluster(cl)

# Check fit
for(i in 1:3) bayesP.ewil_occ <- mean(ewil_occ[, "fit.new",][[i]] > ewil_occ[, "fit",][[i]]) # 0.258 oaky
print(bayesP.ewil_occ, dig = 3)

par(mfrow=c(1,1))
plot(as.matrix(ewil_occ[, "fit",]), as.matrix(ewil_occ[, "fit.new",])) # imperfect but acceptable
abline(0, 1, col = 'red')

print(gelman.diag(x=ewil_occ[,c("alpha.lam", "beta1.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "fit", "fit.new")]), dig=3) # 

gelman.diag(x=ewil_occ)


# convert mcmc list into dataframe
ewil.iter <- as.data.frame(as.matrix(ewil_occ))

# filter to elevation columns
ewil.N.elev <- ewil.iter %>%
  dplyr::select(Nintercept = alpha.lam, Elevation = beta1.lam, Elevation2 = beta2.lam) 

# calculate abundance at each iteration at 
ewil.peak.psi <- rep(NA_real_, times = nrow(ewil.iter))
df.el <- data.frame(el)
for(i in 1:nrow(ewil.iter)) {
  df.el[ , paste0("psi",i)] <- exp(ewil.N.elev[i, "Nintercept"] + ewil.N.elev[i, "Elevation"] * el.s + ewil.N.elev[i, "Elevation2"] * el2.s) / (1 + exp(ewil.N.elev[i, "Nintercept"] + ewil.N.elev[i, "Elevation"] * el.s + ewil.N.elev[i, "Elevation2"] * el2.s))
  ewil.peak.psi[i] <- df.el[which(df.el[ , paste0("psi",i)] == max(df.el[ , paste0("psi",i)])), "el"]
}

# calculate expected peak elevation with 95% CRI
median(ewil.peak.psi)
quantile(ewil.peak.psi, probs = c(0.025, 0.975))

# summarize iterations for plotting
df.el <- df.el %>%
  dplyr::select(-el)
df.el.summary <- data.frame(el)
df.el.summary$median <- apply(df.el, MARGIN = 1, FUN = median)
df.el.summary$LCRI <- apply(df.el, MARGIN = 1, FUN = quantile, probs = c(0.025))
df.el.summary$UCRI <- apply(df.el, MARGIN = 1, FUN = quantile, probs = c(0.975))

# plot elevation vs. expected abundance with CRI
ewil.elev.psi.plot <- ggplot(df.el.summary, aes(el, median)) + geom_line() + xlab("Elevation (m)") + ylab("Probability of occupancy") + ggtitle(expression(italic("E. wilderae"))) + geom_ribbon(aes(ymin=LCRI, ymax=UCRI),alpha=0.3) + theme_bw()

#######################################
# Plethodon jordani - Bayesian Analysis
#######################################
#----------Overdispersion in abundance and detection-----no convergence--------
# random effect of site on abundance and random effect of transect*observation in detection
# Define model
sink("pjor_od2.txt")
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta10.lam*gcover[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]]
    
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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od2.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od3.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

# Results
pjor_od3 <- mcmc.list(out)
plot(pjor_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")]) # 
par(mfrow = c(1,1))
summary(pjor_od3[,c("alpha.lam", "beta1.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta11.lam", "beta13.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p", "beta5.p", "beta10.p", "sigma.site", "sigma.p", "fit", "fit.new")])

stopCluster(cl)

save(pjor_od3, file = "Results/JAGS/pjor_mcmc_out.RData")

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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od4.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od4.txt", pjor.od.data, inits, n.adapt=300000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
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
#----------Overdispersion in abundance and detection reduced model --- no run---------------
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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od5.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od5.txt", pjor.od.data, inits, n.adapt=100, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=1000, thin=1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od2.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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





#----------Overdispersion in detection simple spatial model no elevation^2---not reasonable theta----------

# Define model
sink("pjor_od6.txt")
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta8.lam*ltwi[i] + beta9.lam*canopy[i] + W[i] # + eps.lam[site[i]]
    
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
                     d = d,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od6.txt", pjor.od.data, inits, n.adapt=500, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=1000, thin=1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od2.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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





#----------Overdispersion in detection Gaussian spatial model no elevation^2-------------

# Define model
sink("pjor_od7.txt")
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
    H[i, j] <- (1/tauw) * exp(pow(-d[i,j], 2) / pow(theta, 2)) # Spatial covariance matrix
    }
    muW[i] <- 0
    }
    Omega[1:n.transects, 1:n.transects] <- inverse(H[1:n.transects, 1:n.transects]) # spatial precision mat
    
    # Likelihood
    for(i in 1:n.transects){
    N[i] ~ dpois(lambda[i])
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta8.lam*ltwi[i] + beta9.lam*canopy[i] + W[i] # + eps.lam[site[i]]
    
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
       theta = runif(1, 1, 5),
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
                     d = d,
                     Precip.s = Precip.s[ ,1:6])

library(parallel)
library(rjags)

cl <- makeCluster(4)                       # Request # cores
clusterExport(cl, c("pjor.od.data", "inits", "params", "Nst")) # Make these available
clusterSetRNGStream(cl = cl, 1259)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od6.txt", pjor.od.data, inits, n.adapt=500, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=1000, thin=1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("pjor_od2.txt", pjor.od.data, inits, n.adapt=100000, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=100000, thin=100) # Sample from posterior distribution
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
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
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
    out <- coda.samples(jm, params, n.iter=100000, thin=40) # Sample from posterior distribution
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


###################################################
# D. ocoee & imitator combined - Bayesian Analysis
###################################################

Dimioco <- DIMI + DOCO

#---------------Dimioco Reduced overdispersion model-no convergence-------------
# Define model
sink("dimioco_od1.txt")
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
    #beta10.lam ~ dnorm(0, 0.01)
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta2.lam*elev2[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta10.lam*gcover[i]
    
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

Nst <- apply(Dimioco, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       beta2.lam = rnorm(1, -1, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0.2, 1),#,
       sigma.p = runif(1, 0.2, 1))#,
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
             #"beta10.lam",
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
             #"delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")#,
#"eval",
#"y.new")

dimioco.od.data <- list(C = as.matrix(Dimioco[, 1:6]), 
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
                     Temp.s = Temp.s[ ,1:6],
                     Temp.s2 = Temp.s*Temp.s,
                     RH.s = RH.s[ ,1:6],
                     Precip.s = Precip.s[ ,1:6])

runif(1)

CL <- makeCluster(4) # set nc (usually 3) cores
clusterExport(cl = CL, list("dimioco.od.data", "params", "inits", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)

system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dimioco_od1.txt", dimioco.od.data, inits, n.adapt = 50000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dimioco.od1.list <- mcmc.list(out)
plot(dimioco.od1.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")]) # seems good, maybe slight autocorr
par(mfrow = c(1,1)) # one bad chain

plot(out[[1]][,"sigma.p"]) # bad chain
plot(out[[2]][,"sigma.p"])
plot(out[[3]][,"sigma.p"])
plot(out[[4]][,"sigma.p"])

dimioco.od1.list <- mcmc.list(out[[3]], out[[2]], out[[4]])
summary(dimioco.od1.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")])
plot(dimioco.od1.list[,c("alpha.lam", "beta1.lam", "beta2.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")]) # seems good, maybe slight autocorr
par(mfrow = c(1,1)) # no convergence of elev and elev2

stopCluster(CL)


#---------------Dimioco Reduced overdispersion model-without elev2-------------
# Define model
sink("dimioco_od2.txt")
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
    #beta10.lam ~ dnorm(0, 0.01)
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
    
    log(lambda[i]) <- alpha.lam + beta1.lam*elev[i] + beta3.lam*slope[i] + beta4.lam*aspectN[i] + beta5.lam*aspectE[i] + beta6.lam*tpi[i] + beta7.lam*trail[i] + beta8.lam*ltwi[i]+ beta9.lam*canopy[i] + beta11.lam*litterdepth[i] + beta12.lam*slope2[i] + beta13.lam*lstream[i] + eps.lam[site[i]] #+ beta10.lam*gcover[i]
    
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

Nst <- apply(Dimioco, 1, function(x) max(x, na.rm = TRUE)) + 1
inits <- function(){
  list(N = Nst,
       alpha.lam = rnorm(1, -1, 1),
       beta1.lam = rnorm(1, 2.5, 1),
       #beta2.lam = rnorm(1, -1, 1),
       beta3.lam = rnorm(1, 0, 1),
       beta4.lam = rnorm(1, 0, 1),
       beta5.lam = rnorm(1, 0, 1),
       beta6.lam = rnorm(1, 0, 1),
       beta7.lam = rnorm(1, 0, 1),
       beta8.lam = rnorm(1, 0, 1),
       beta9.lam = rnorm(1, 0, 1),
       # beta10.lam = rnorm(1, 0, 1),
       beta11.lam = rnorm(1, 0, 1),
       beta12.lam = rnorm(1, 0, 1),
       beta13.lam = rnorm(1, 0, 1),
       eps.lam = rnorm(48, 0, 1),
       alpha.p = rnorm(1, -1, 1),
       beta1.p = rnorm(1, 0, 1),
       beta2.p = rnorm(1, 0, 1),
       beta3.p = rnorm(1, -0.5, 1),
       beta4.p = rnorm(1, 0.5, 1),
       beta5.p = rnorm(1, 0.5, 1),
       beta10.p = rnorm(1, 0.5, 1),
       sigma.site = runif(1, 0.2, 1),#,
       sigma.p = runif(1, 0.2, 1))#,
  #delta.p = rnorm(195*6, 0, 1))
}

params <- c( "alpha.lam", 
             "beta1.lam", 
            # "beta2.lam", 
             "beta3.lam",
             "beta4.lam",
             "beta5.lam",
             "beta6.lam",
             "beta7.lam",
             "beta8.lam",
             "beta9.lam",
             #"beta10.lam",
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
             #"delta.p",
             "sigma.site",
             "sigma.p",
             "N",
             #"p",
             "fit",
             "fit.new")#,
#"eval",
#"y.new")

dimioco.od.data2 <- list(C = as.matrix(Dimioco[, 1:6]), 
                        n.transects = n.transects, 
                        n.surveys = n.surveys,
                        n.sites = length(site.inits), 
                        elev = elev, 
                       # elev2 = elev2, 
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

runif(1)

CL <- makeCluster(6) # set nc (usually 3) cores
clusterExport(cl = CL, list("dimioco.od.data2", "params", "inits", "Nst")) # make data available to each core
clusterSetRNGStream(cl = CL, 4345)

system.time({
  out <- clusterEvalQ(CL, {
    library(rjags)
    jm <- jags.model("dimioco_od2.txt", dimioco.od.data2, inits, n.adapt = 50000, n.chains = 1)
    fm <- coda.samples(jm, params, n.iter = 5000, thin = 5)
    return(as.mcmc(fm))
  }) 
}) # 

# Results
dimioco.od2.list <- mcmc.list(out)
stopCluster(CL)

plot(dimioco.od2.list[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "alpha.p", "beta1.p", "beta2.p", "sigma.site", "sigma.p")]) # on bad chain
par(mfrow = c(1,1)) # 

plot(out[[1]][,"alpha.p"]) # 
plot(out[[2]][,"alpha.p"]) # bad chain
plot(out[[3]][,"alpha.p"])
plot(out[[4]][,"alpha.p"])

dimioco.od2.list <- mcmc.list(out[[1]], out[[3]], out[[4]], out[[5]], out[[6]])
summary(dimioco.od2.list[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta13.lam", "beta11.lam", "beta12.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p",  "beta4.p",  "beta5.p",  "beta10.p",  "sigma.site", "sigma.p")])
plot(dimioco.od2.list[,c("alpha.lam", "beta1.lam", "beta3.lam", "beta4.lam", "beta5.lam", "beta6.lam", "beta7.lam", "beta8.lam", "beta9.lam", "beta13.lam", "beta11.lam", "beta12.lam", "alpha.p", "beta1.p", "beta2.p", "beta3.p",  "beta4.p",  "beta5.p",  "beta10.p",  "sigma.site", "sigma.p")]) # acceptable
par(mfrow = c(1,1)) # 

stopCluster(CL)


N.dimioco.median <- summary(dimioco.od2.list[,c("N[1]")])$quantiles["50%"]
N.dimioco.median$quantiles["50%"]

N.dimioco.mean <- matrix(NA, 195, 1)
N.dimioco.uci <- matrix(NA, 195, 1)
N.dimioco.lci <- matrix(NA, 195, 1)
for(i in 1:195){
  foo[i] <- paste("N[",i,"]", sep="")
  bar <- summary(dimioco.od2.list[,c(foo[i])])$statistics["Mean"]
  attributes(bar) <- NULL
  N.dimioco.mean[i] <- bar
  N.dimioco.uci[i] <- summary(dimioco.od2.list[,c(foo[i])])$quantiles["97.5%"]
  N.dimioco.lci[i] <- summary(dimioco.od2.list[,c(foo[i])])$quantiles["2.5%"]
}

DIMIOCOmin <- apply(Dimioco, 1, function(x) max(x, na.rm = TRUE))

png("OcoeeAbundance.png", width=6, height=6, units='in', res = 150)
plot(Data$elev, N.dimioco.mean, pch=16, ylab = "Abundance of occee complex", xlab = 'Elevation (m)')
points(Data$elev, DIMIOCOmin, col = 'red') 
legend(400, 90, legend=c("Abundance", "Max Count"), pch = c(16, 1), col = c(1, 2))# bimodal likely separates imitator and ocoee but could also just be due to low captures of both species in shitty mid-elevation rhododendron sites
dev.off()




