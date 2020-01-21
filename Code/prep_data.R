#########################################################################
# GSMNP Elevation Gradient Analysis
# October 31, 2012
# Daniel J. Hocking
########################################################################

#-------Load Packages------
# library(parallel)
# library(rjags)
library(dplyr)
library(ggplot2)
library(devtools)
# install_github(repo = "ropensci/prism")
library(prism) ##prism data access
library(tidyr)
library(stringr)

#-----------Load Data and do Manipulation--------
# Import and check data
Data <- read.table('Data/Raw/TransectData2012.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
Side <- read.table('Data/Raw/Side.csv', header = TRUE, sep = ',')
Stream <- read.table('Data/Raw/Strm_Dist.csv', header = TRUE, sep = ",")
site.inits <- read.table('Data/Raw/site_inits.csv', header = TRUE, sep = ',')

# merge data
Data <- merge(x=Data, y=Side, by.x = "site", by.y = 'site', all.x = TRUE, sort = FALSE)
Data <- merge(x = Data, y=Stream, by.x = "transect", by.y = 'plot_trans', all.x = TRUE, sort = FALSE)

# get number of transects
n.transects <- nrow(Data)

# Make numeric and factor versions of site groupings
# Data <- Data %>%
#   dplyr::arrange(site)
# 
# Data$site1 <- as.numeric(seq(1,n.transects))
# for(i in 2:n.transects){
#   Data$site1[i] <- ifelse(Data$site[i] == Data$site[i-1], Data$site1[i-1], (Data$site1[i-1] + 1))
# }

# Data$sitef <- as.factor(Data$site1)
Data$sitef <- as.factor(Data$site)
sites <- unique(Data$site)
sites_df <- data.frame(site = sites, site_num = seq(1, length(sites)), stringsAsFactors = FALSE)
Data <- Data %>%
  dplyr::left_join(sites_df)
Data$sidef <- as.factor(Data$side)

# Data <- Data[order(Data$elev, na.last = TRUE), ] 

# separate data
PJOR <- Data[ , c('PJOR1', 'PJOR2', 'PJOR3', 'PJOR4', 'PJOR5', 'PJOR6')]
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

PJOR <- replace(PJOR, is.na(RH), NA)
DCON <- replace(DCON, is.na(RH), NA)
DIMI <- replace(DIMI, is.na(RH), NA)
DMON <- replace(DMON, is.na(RH), NA)
DOCO <- replace(DOCO, is.na(RH), NA)
DQUA <- replace(DQUA, is.na(RH), NA)
DSAN <- replace(DSAN, is.na(RH), NA)
DWRI <- replace(DWRI, is.na(RH), NA)
EWIL <- replace(EWIL, is.na(RH), NA)
GPOR <- replace(GPOR, is.na(RH), NA)
NVIR <- replace(NVIR, is.na(RH), NA)
PGLU <- replace(PGLU, is.na(RH), NA)
PSER <- replace(PSER, is.na(RH), NA)
PTEY <- replace(PTEY, is.na(RH), NA)

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

write.table(data.frame(Species, Count), file = 'Results/TotalCounts.csv', sep = ",", row.names = FALSE)

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

#-------------- Explore Data ---------------------
# Naive examination without accounting for detection
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

#---------------- Prep Data for Full Analysis -----------------
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
# Temp2.s <- (Temp2 - mean(Temp2, na.rm = TRUE))/sd(Temp2, na.rm = TRUE)

# PLJO <- as.matrix(PLJO)
# dimnames(PLJO) <- NULL

Rain <- as.matrix(Rain + 1)
dimnames(Rain) <- NULL

VegWet <- as.matrix(VegWet + 1)
dimnames(VegWet) <- NULL

Ground <- as.matrix(Ground + 1)
dimnames(Ground) <- NULL

tpi <- scale(as.vector(Data$tpi_100)) # negative lower than surrounding (valley), positive (ridge), 0 is flat or constant slope
attributes(tpi) <- NULL
str(tpi)

solar <- (Data$solar - mean(Data$solar))/sd(Data$solar)
twi <- (Data$twi_10 - mean(Data$twi_10))/sd(Data$twi_10)
elev <- (Data$elev - mean(Data$elev))/sd(Data$elev)
elev2 <- elev^2
slope <- (Data$Slope- mean(Data$Slope))/sd(Data$Slope)
slope2 <- slope^2
Data$AspectR <- Data$Aspect*2*pi/360 # convert degrees to radians
Data$aspectN <- cos(Data$Aspect)
Data$aspectE <- sin(Data$Aspect)
aspectN <- Data$aspectN
aspectE <- Data$aspectE
trail <- Data$trail
site <- Data$sitef
canopy <- (Data$canopy - mean(Data$canopy))/sd(Data$canopy)
litterdepth <- (Data$litdepth - mean(Data$litdepth))/sd(Data$litdepth)
gcover <- (Data$gcover - mean(Data$gcover, na.rm = TRUE))/sd(Data$gcover, na.rm = TRUE)
Data$ltwi <- log(Data$twi_10 + 10)
# ltwi <- Data$twi_10
ltwi <- (Data$ltwi - mean(Data$ltwi))/sd(Data$ltwi)
Data$lstream <- log(Data$strm_dist + 1)
lstream <- (Data$lstream - mean(Data$lstream))/sd(Data$lstream)
stream <- (Data$strm_dist - mean(Data$strm_dist))/sd(Data$strm_dist)
gcover2 <- gcover^2

site.inits <- as.vector(site.inits$siteinit)

n.surveys <- ncol(PJOR)

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

d_m <- as.matrix(dist(sp::coordinates(res@coords[ , 1:2]))) # columns 1:2 are the x,y coordinates (UTM) - only problem with this method is keeping transect names. Will work if order doesn't change anywhere.

d_km <- d_m / 1000 # convert distances to kilometers


#---------------- Prism Precipitation ---------------------
if(!dir.exists("Data/Prism/Precip")) dir.create("Data/Prism/Precip", recursive = TRUE)
options(prism.path = "Data/Prism/Precip")
get_prism_normals(type = 'ppt', resolution = '800m', mon = 5:10, keepZip = TRUE)
mystack <- ls_prism_data(absPath = FALSE) %>%  prism_stack()  # Stack files
 
# Get proj from raster stack
mycrs <- mystack@crs@projargs
#original_crs <- xy@proj4string

library(sp)

# Convert points to spatial points data frame
# coordinates(xy) <- c('lon', 'lat')
xy1 <- spTransform(xy, CRS(mycrs))
proj4string(xy1) <- CRS(mycrs)

xy <- Data[ , c("Long", "Lat")]
coordinates(xy) <- c('Long', 'Lat')
proj4string(xy1) <- CRS(mycrs)
xy1 <- spTransform(xy, CRS(mycrs))

# Extract data from raster
ppt_30yr <- data.frame(coordinates(xy1), transect = Data$transect, raster::extract(mystack, xy1), stringsAsFactors = FALSE)

# Reshape data
ppt_30yr <- ppt_30yr %>%  gather(month, value, 4:ncol(ppt_30yr))

# Split date into type, year, month, and day
ppt_30yr$month <- gsub("PRISM_ppt_30yr_normal_800mM2_", "", ppt_30yr$month) %>% gsub("_bil", "", .) %>% as.numeric() 

# Reshape data again
# ppt_30yr <- ppt_30yr %>%  spread(type, value)

ppt_30yr <- dplyr::rename(ppt_30yr, ppt = value)
str(ppt_30yr)

# Order data
# ppt_30yr <- ppt_30yr[order(ppt_30yr$transect),]

str(ppt_30yr)
head(ppt_30yr)

# summarize
ppt_30yr_active <- ppt_30yr %>%
  # dplyr::rename(lon = Long, lat = Lat) %>%
  dplyr::group_by(lon, lat, transect) %>%
  dplyr::summarise(ppt = min(ppt)) # use month with minimum ppt at each transect

str(ppt_30yr_active)
summary(ppt_30yr_active)

# plot(ppt_30yr_active)

plot(mystack[[1]], xlim = as.numeric(xy1@bbox[1, ]), ylim = as.numeric(xy1@bbox[2, ]), main = "Monthly Precip")
points(xy1@coords)
# points(xy1@coords) # add color to represent elevation? Add hillshade or contour lines?

#---------------- Prism Temperature ---------------------

if(!dir.exists("Data/Prism/Temp")) dir.create("Data/Prism/Temp", recursive = TRUE)

options(prism.path = "Data/Prism/Temp")
get_prism_normals(type = 'tmean', resolution = '800m', mon = 5:10, keepZip = TRUE)
temp_stack <- prism_stack(ls_prism_data(absPath = FALSE)) # Stack files

# Get proj from raster stack
mycrs <- temp_stack@crs@projargs
#original_crs <- xy@proj4string

# Convert points to spatial points data frame
# coordinates(xy) <- c('lon', 'lat')
xy <- Data[ , c("Long", "Lat")]
coordinates(xy) <- c('Long', 'Lat')
proj4string(xy1) <- CRS(mycrs)
# xy1 <- spTransform(xy, CRS(mycrs))

# Extract data from raster
temp_30yr <- data.frame(coordinates(xy1), transect = Data$transect, raster::extract(temp_stack, xy1), stringsAsFactors = FALSE) # transect = dimnames(xy@coords)[[1]]

# Reshape data
temp_30yr <- temp_30yr %>%  gather(month, value, 4:ncol(temp_30yr))

# Split date into type, year, month, and day
temp_30yr$month <- gsub("PRISM_tmean_30yr_normal_800mM2_", "", temp_30yr$month) %>% gsub("_bil", "", .) %>% as.numeric() 

# Reshape data again
# temp_30yr <- temp_30yr %>%  spread(type, value)

temp_30yr <- dplyr::rename(temp_30yr, tmean = value)
str(temp_30yr)

# Order data
# temp_30yr <- temp_30yr[order(temp_30yr$transect),]

str(temp_30yr)
head(temp_30yr)

# summarize
temp_30yr_active <- temp_30yr %>%
  # dplyr::rename(lon = Long, lat = Lat) %>%
  dplyr::group_by(lon, lat, transect) %>%
  dplyr::summarise(tmean = max(tmean)) # maybe use hottest month? rather than mean

str(temp_30yr_active)
summary(temp_30yr_active)

# plot(temp_30yr)

# standardize
temp_30yr_active$temp30 <- (temp_30yr_active$tmean - mean(temp_30yr_active$tmean)) / sd(temp_30yr_active$tmean)
ppt_30yr_active$ppt30 <- (ppt_30yr_active$ppt - mean(ppt_30yr_active$ppt)) / sd(ppt_30yr_active$ppt)


# Combine and summarize data
df <- Data %>%
  select(transect, site, Lat, Long, tpi_100, solar, twi_10, Slope, aspectN, aspectE, canopy, trail, elev) %>%
  left_join(temp_30yr_active, by = "transect") %>%
  left_join(dplyr::select(ungroup(ppt_30yr_active), -lat, -lon), by = "transect")
summary(df)

library(ggplot2)
df$elev_bin <- as.numeric(cut_number(df$elev, 5))

plot(temp_stack[[1]], xlim = as.numeric(xy1@bbox[1, ]), ylim = as.numeric(xy1@bbox[2, ]), main = "Monthly Temperature")
points(xy1@coords, pch = 16, col = df$elev_bin)
legend("topright", "Pt color = Elevation")
# points(xy1@coords) # add color to represent elevation? Add hillshade or contour lines?


#--------------- Check for correlations ------------------
# Scatterplot matrix
source('/Users/Dan/Documents/Statistics/R/Functions/panelcor.R') # get from github utility functions eventually
Pairs1 <- data.frame(DWRImin, PJORmin, Data$elev, Data$tpi_100, Data$trail, Data$aspectN, Data$aspectE, Data$Slope, Data$solar, Data$twi_10)

pairs(Pairs1, upper.panel=panel.smooth, lower.panel=panel.cor) # nothing excessively correlated

Pairs2 <- data.frame(PJORmin, temp_30yr_active$tmean, ppt_30yr_active$ppt, Data$elev, Data$tpi_100, Data$Slope, Data$solar, Data$twi_10, canopy)

pairs(Pairs2, upper.panel=panel.smooth, lower.panel=panel.cor) # something not right!!!!!

pairs(select(df, -transect, -site, -Lat, -Long, -elev_bin, -lon, -lat, -trail), upper.panel=panel.smooth, lower.panel=panel.cor)

#----------- remove NA in independent variables -------------

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)) # should be a problem since no data in these dependent variables either
Temp.s <- impute.mean(Temp.s)
RH.s <- impute.mean(RH.s)
Precip.s <- impute.mean(Precip.s)

gcover <- impute.mean(gcover)
gcover <- impute.mean(gcover)
gcover2 <- gcover^2

gcover2 <- gcover^2
Temp.s2 <- Temp.s^2

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

# Summary
summary(Data)

#------------- organize and summarize transect info --------------
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

#---------------- summarize elevational info -------------
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
write.table(min.max.elev, file = "Results/min-max-elev.csv", sep = ",")


#------------ save everything for use in models -----------

if(!dir.exists("Data/Derived")) dir.create("Data/Derived", recursive = TRUE)

n.sites <- length(unique(Data$site))
n.sites == max(Data$site_num)

save(Data, df, Precip.s, Temp.s, RH.s, solar, twi, elev, slope, aspectE, aspectN, trail, site, canopy, gcover, litterdepth, ltwi, stream, d_m, d_km, n.transects, n.sites, n.surveys, PJOR, PJORmin, DWRI, DWRImin, EWIL, EWILmin,  file = "Data/Derived/jags_prep.RData")



## Remove missing observations for now
summary(PJOR)
PJOR5_na <- PJOR[ , 1:5]

na_rows <- which(is.na(rowSums(PJOR5_na)))
rows <- which(!is.na(rowSums(PJOR5_na)))
PJOR5 <- PJOR5_na[which(!is.na(rowSums(PJOR5_na))), ]
dim(PJOR5)
summary(PJOR5)

elev5 <- elev[rows]
stream5 <- stream[rows]
twi5 <- twi[rows]
litter5 <- litterdepth[rows]

RH5 <- RH.s[rows, 1:5]
temp5 <- Temp.s[rows, 1:5]
precip5 <- Precip.s[rows, 1:5]
gcover5 <- gcover[rows]

n.transects <- length(elev5)
n.surveys <- ncol(PJOR5)

Data5 <- Data[rows, ]
Data5 <- Data5 %>%
  arrange(site, transect)
Data5$site_stan <- 99999
Data5$site_stan[1] <- 1
for(i in 2:n.transects) {
  Data5$site_stan[i] <- ifelse(Data5$site[i] == Data5$site[i-1], Data5$site_stan[i-1], Data5$site_stan[i-1] + 1)
}

n.sites <- length(unique(Data5$site_stan))

save(Data, Data5, Precip.s, Temp.s, RH.s, solar, twi, elev, slope, aspectE, aspectN, trail, site, canopy, gcover, litterdepth, ltwi, stream, d_m, d_km, n.transects, n.sites, n.surveys, PJOR, PJOR5, elev5, stream5, twi5, litter5, RH5, temp5, precip5, gcover5, n.transects, n.surveys, n.sites,  file = "Data/Derived/stan_prep.RData")



save.image("Data/Derived/prep_image.RData")

rm(list = ls())
