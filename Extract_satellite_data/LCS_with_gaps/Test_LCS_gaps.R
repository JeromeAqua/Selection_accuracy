args <- commandArgs(TRUE)
# J <- strtoi(args[1], base=10L)


jj <- strtoi(args[1], base=10L)#floor(J[1]/4) + 1 # This is the variable we work with -- 37 for sst
# 
# tt <- (J[1] %% 4) 


library("ncdf4")
library("tidyverse")
library("glue")
library("lubridate")
# library("pracma") #for interp2 and meshgrid
# setwd("/home/jpinti/Extract_satellite_data/")
setwd("/home/2453/Importance_accuracy")

source("../Extract_satellite_data/interp2.R") #to returnv NA when there is NA in query points -- not the case with the pracma package
source("../Extract_satellite_data/meshgrid.R")

library("data.table") #for fcoalesce
# library("sendmailR")

# library("foreach")
# library("doParallel")
library("parallel")
## Load tracking data
# load("all_tracks0B.RData")
# Trackingdates <- unique(Synthetic$date)
# Trackingdates <- ymd_hms(Trackingdates)
# Trackingdates <- Trackingdates[as.numeric(format(Trackingdates, format = "%Y")) >2001] ## Remove data points before 2001 as we won't be able to interpolate anyways

# Limits$XMIN[LIMITS$XMIN < -180] <- -180
# Limits$XMAX[LIMITS$XMAX > 180] <- 180

## Load satellite data
LCS <- nc_open("/lustre/scratch/orb/lsq/daily/lsq_daily_8km_2006.nc4", write = FALSE)
ModisAquaDaily <- nc_open("/lustre/scratch/orb/aqua/8day/aqua_8day.2006.concat.nc4", write = FALSE)
yr <- 2006
#"http://basin.ceoe.udel.edu/thredds/dodsC/AQUAGLOBAL8DAY9KM.nc")
# setwd("/home/2453/Extract_satellite_data/")

Latitude <- ncvar_get(LCS, varid = "latitude", start = 1, count = -1)
Longitude <- ncvar_get(LCS, varid = "longitude", start = 1, count = -1)
Time <- ncvar_get(LCS, varid = "time", start = 1, count = -1)
Time <- as.POSIXct(Time, origin = "1970-01-01", tz = "UTC")

LatitudeCHL <- ncvar_get(ModisAquaDaily, varid = "lat", start = 1, count = -1)
LongitudeCHL <- ncvar_get(ModisAquaDaily, varid = "lon", start = 1, count = -1)
TimeCHL <- ncvar_get(ModisAquaDaily, varid = "time", start = 1, count = -1)
TimeCHL <- as.POSIXct(TimeCHL, origin = "1970-01-01", tz = "UTC")


ftle <- ncvar_get(LCS, varid = "ftle", start = c(1,1,30), count = c(-1,-1,1))  # then we extract what we need from the satellites
chl <- ncvar_get(ModisAquaDaily, varid = "chl_ocx", start = c(1,1,30), count = c(-1,-1,1))  # then we extract what we need from the satellites


Products <- names(LCS$var)
nlon <- length(Longitude)

Products
# print(glue("{Products[jj]}"))

nc_close(LCS)
nc_close(ModisAquaDaily)
warnings()

# #stop cluster use when we are done
# parallel::stopCluster(cl = my.cluster)

save.image("test_lcs_gaps.RData")
# save(outT, outB, outC, outL, outJ, file = paste(glue("/lustre/scratch/jpinti/Accuracy/Accuracy_lcs_matched_{jj}_111.RData")))

####################################
MGG <- meshgrid(Longitude, Latitude)
MGGlongs <- c(MGG$X)
MGGlats <- c(MGG$Y)
library(tictoc)
tic()
chl2 <- interp2(LatitudeCHL, LongitudeCHL, chl, MGGlats, MGGlongs, method = "nearest")
toc()

chl3 <- matrix(chl2, nrow = length(Longitude), byrow = TRUE)

ftle[is.na(chl3)] <- NA
