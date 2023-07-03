library(sp)
require(maptools)
#library(tictoc)

library("tidyverse")
library("ncdf4")
library("maps")
library("sf")
library("raster")
library("rgdal")

library("dplyr")
library("geosphere")

#library("marmap")
#library("spData")
#library("spDataLarge")
#library("tmap")    # for static and interactive maps
#library("leaflet") # for interactive maps
library("circular")

library("ggplot2") # tidyverse data visualization package
library("sendmailR")

library("foreach")
library("doParallel")

library("lubridate")


# library(ggnewscale)
# library(RColorBrewer)

# Function to know if a point is in the ocean
point.in.SpatialPolygons = function(point.x, point.y, SpP){
  # VERBATIM COPY FROM prevR package. They deserve all the credit for this function
  ###############################################################################################
  # Cette fonction renvoie pour chaque point defini par le couple (point.x, point.y) T ou F 
  #     si le point est a l'interieur ou non du spatialPolygons SpP
  # Un point est considere a l'interieur de N polygons si il est a l'interieur d'au moins
  # un polygon non Hole et a l'exterieur de tous les polygons Hole
  # Cette foncion est utilisee par toutes les fonctions de lissage (krige , kde, idw) . 
  # En effet ces fonctions travaillent sur un grid rectangulaire englobant les donnees. En presentation on ne veut que les resulats 
  #   interieurs a la frontiere qui est definie dans l'element SpP (SpP contient boundary). 
  #   Tous les elements du grid hors de la frontiere seront dans les programmes de lissage positonnes a NA
  # 
  ###############################################################################################
  
  X = slot(SpP,"polygons")
  is.inside = F
  for(i in 1:length(X)){
    PS   = slot(X[[i]],"Polygons")
    for(j in 1:length(PS)){
      pol.x = slot(PS[[j]],"coords")[,1]
      pol.y = slot(PS[[j]],"coords")[,2]
      pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
      if(!slot(PS[[j]],"hole")) {
        is.inside = is.inside | pointsPosition != 0
      }
    }
  }
  is.outsideHole = T
  for(i in 1:length(X)){
    PS   = slot(X[[i]],"Polygons")
    for(j in 1:length(PS)){
      pol.x = slot(PS[[j]],"coords")[,1]
      pol.y = slot(PS[[j]],"coords")[,2]
      pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
      if(slot(PS[[j]],"hole")) {
        is.outsideHole = is.outsideHole & (pointsPosition == 0 |  pointsPosition == 3)
      }
    }
  }
  is.inside & is.outsideHole
}



# Read the oceans' shapefile
setwd("/home/jpinti/Importance_accuracy/") # C:/Users/jpinti/Documents/R/Ocean_mask/
# load("../Reanalysis2.RData")
oceans <- readOGR("ne_10m_ocean.shp")
# xcarib <- c(-104.6316, -95.5350, -89.91, -84.5047, -81.9189, -79.1083, -77.0419, -68.48, -60, -104.6316)
# ycarib <- c(25.9368, 17.0343, 16.3186, 11.3277, 8.551979, 9.229575, 7.6355, 8.8436, 28, 25.9368)
# Carib1 <- Polygon(cbind(xcarib, ycarib))
# Carib2 <- Polygons(list(Carib1), ID = "A")
# Carib <- SpatialPolygons(list(Carib2))

# issues <- as.character(scan("Potentialteleportingleatherbacks.txt"))
# issues <- c(issues, unique(TOPP$event_id[TOPP$species_id==70]), 720600500, 1410)
# issues <- unique(c(issues, 220600400, 1044, 601, 700701300, 700500300, 700500900))

# Initialize list of datasets
# events <- unique(TOPP$event_id)
# TOPP$id <- seq.int(nrow(TOPP)) #1 number for each line of the dataset

# Create the cluster for parallel computing
# my.cluster <- parallel::makeCluster(
#   50,
#   type = "PSOCK")
# doParallel::registerDoParallel(cl = my.cluster)
# # cl <- makeCluster(30)
# # registerDoParallel(cl)  # use multicore, set to the number of our cores
# 
# 
# ################################################################################################################
# ################################ Start RANDOM WALK HERE ######################################################
# ################################################################################################################
# 
# 
# load("Syntheticpart_nooceantest.RData")
nc <- nc_open("MODIS_8day_chl_ocx.nc") #"http://basin.ceoe.udel.edu/thredds/dodsC/AQUAGLOBAL8DAY9KM.nc") #"MODIS_8day_sst.nc")

# Synthetic = foreach (runnumber=1:100, .export = "point.in.SpatialPolygons", .packages =
#                 c("tidyverse","circular","dplyr","geosphere","sp","lubridate","ncdf4"), .combine = "c" ) %dopar% {

TOT <- list()
print("starting")
for (runnumber in 1:100){
  
  
  # concentration <- 5
  set.seed(2*runnumber) 
  for (concentration in 0){#c(20,10,5,2,1,0.75,0.5,0.25,0)){########################################################################
    
    #initialization
    StartD <- as.POSIXct(strptime("2006-05-01 12:00:00", "%Y-%m-%d %H:%M:%S"), origin = "1970-01-01", tz = "UTC")
    Startlat <- runif(1)*30+10 #between 10 and 40
    Startlon <- runif(1)*10-140 #between -140 and -130
    durtrack <- 80
    
    a <-0
    Simulationstart <- data.frame(event_id = runnumber, kappa = concentration, date = as.character(ymd_hms(StartD)), dt = 1, latitude = Startlat, longitude = Startlon, bearing = runif(1)*360, steplength = NA)
    
    if (concentration==0){Simulation <- Simulationstart
    Simulation$date <- as.character(Simulationstart$date)} else {Simulation[nrow(Simulation)+1,] <- Simulationstart
    Simulation$date <- as.character(Simulationstart$date)}
    
    newp <- c(Startlon, Startlat)
    
    # #For the distance if non biased walk
    STEP <- function(kappa) {1000*as.vector(rnorm(1,0,10/(kappa+0.5)))} #in meters
    #For the angle
    ANGLE <- function(kappa) {as.vector(rvonmises(1, circular(0), kappa, control.circular=list(units="degrees")))} #kappa is the concentration, the higher it is the more the track is straight
    
    ##For the sst bias
    lon <- ncvar_get(nc, varid = "longitude")
    lat <- ncvar_get(nc, varid = "latitude")
    Time <- ncvar_get(nc, varid = "time", start = 1, count = -1)
    Time <- as.POSIXct(Time, origin = "1970-01-01", tz = "UTC")
    
    for(i in 2:durtrack){ 
      # Computation of the move due to SST
      # print(i)
      satellitetime <- which.min(abs(Time - ymd_hms(Simulation$date[nrow(Simulation)]))) #the satellite time we are going to use
      whichlon <- which.min(abs(as.numeric(Simulation$longitude[nrow(Simulation)]) - lon)) #lonlat of the position in the sst matrix
      whichlat <- which.min(abs(as.numeric(Simulation$latitude[nrow(Simulation)]) - lat))
      
      SST <- t(ncvar_get(nc, varid = "chl_ocx", start = c(1,1,satellitetime), count = c(-1,-1,1))) #get the chl ocx of the correct day -still called sst because easier not to change the name
      
      latmini <- (whichlat-5):(whichlat+5) #only look within 50km of the position
      lonmini <- (whichlon-5):(whichlon+5)
      
      SST <- SST[latmini, lonmini] #sst within 50km of the position
      
      if (sum(SST, na.rm = T)>0){   #to avoid bugs when no data
        maxsst <- which(SST == max(SST, na.rm = T),  arr.ind = TRUE)
        maxsstlon <- lonmini[maxsst[2]] #position of the max-min sst within 50km
        maxsstlat <- latmini[maxsst[1]]
        
        pos <- c(as.numeric(Simulation$longitude[nrow(Simulation)]), as.numeric(Simulation$latitude[nrow(Simulation)]))
        maxpos <- c(lon[maxsstlon], lat[maxsstlat])
        
        B <- bearing(pos,maxpos)
        D <- distHaversine(pos,maxpos)
        
        B <- B + ANGLE(concentration) ## Change strength of the selection, in degrees -- % 5 for high selection, 1.5 for low selection
        D <- STEP(concentration)##D + STEP(concentration)#Change strength of the selection, in m -- % 2000 for high selection, 10000 for low selection####################################################
      } else {
        #if all SST is NaN, choose random angle and step length
        B <- ANGLE(concentration) 
        D <- STEP(concentration)
      }
      newp <- destPoint(newp,B,D) #minus angle because destPoint uses clockwise orientation
      
      # oceantest <- point.in.SpatialPolygons(newp[1], newp[2], oceans) #example of test to see if lon / lat is in the ocean
      
      comptbreak <- 0
      
      # while (oceantest==FALSE  & comptbreak < 100){ 
      #   #if the point is not in the ocean, try again for a while
      #   B <- B + rnorm(1, 0, sd = 5) ## Move a bit the point
      # 
      #   newptempo <- destPoint(newp,B,STEP)   
      #   oceantest <- point.in.SpatialPolygons(newptempo[1], newptempo[2], oceans) #example of test to see if lon / lat is in the ocean
      #   comptbreak <- comptbreak + 1
      # }
      # 
      # if (comptbreak > 0){
      #   newp <- newptempo
      #   angle <- 0 #no direction change if unsuccessful time step simulation
      # }
      # 
      # if (comptbreak==100){
      #   # Break for loop if didn't reach the sea
      #   print(c("!!!stuck on land"))
      #   newp <- c(NA, NA)
      #   #so that we can restart a new track at the next time step
      #   # if (Simulation$event_id[i+1]==Simulation$event_id[i]){
      #   #   a <- a-1 #because we will increase it of 1 at the next time step
      #   # }
      #   # save.image(file = "CRW_emp5.RData")
      #   #      stop("stuck on land")
      #   break
      # }
      
      # Then update the position that is in the sea
      newd <- ymd_hms(Simulation$date[nrow(Simulation)])+lubridate::days(1)
      Simulation[nrow(Simulation)+1,] <- c( runnumber, concentration, as.character(newd)  , 1, 
                                            as.numeric(newp[2]), as.numeric(newp[1]), B, D/1000)
      # save.image("tempo.RData")
      # Simulation <- bind_rows(Simulation, TMP)
      
    }
  }
  
  
  # list(Simulation)
  TOT[[runnumber]] <- Simulation
  save.image(file = "Fakeanimals_chl_ocx_0.RData")
  print(runnumber)
}
nc_close(nc)

# # # # #Plot all tracks
# library("marmap")
# library("ggnewscale")
# Simulation$longitude <- as.numeric(Simulation$longitude)
# Simulation$latitude <- as.numeric(Simulation$latitude)
# Simulation$kappa <- as.numeric(Simulation$kappa)
# Pacific<- getNOAA.bathy(lon1=-180,lon2=180,lat1=-66,lat2=67, resolution = 10, keep = TRUE)
# 
# xlims <- range(c(min(Simulation$longitude,na.rm = T), max(Simulation$longitude, na.rm = T)))#c(-65 %%360,130))
# ylims <- range(c(min(Simulation$latitude, na.rm = T), max(Simulation$latitude, na.rm = T)))#c(-50, 65))
# P <-  ggplot(subset(Simulation, Simulation$kappa==0), aes(x = longitude %%360, y = latitude)) + geom_tile(data = Pacific, aes(x = x %%360, y = y, fill = z)) +
#   scale_fill_gradientn(colors = c("darkblue", "lightblue",grey(.7), grey(.9)),values = scales::rescale(c(min(Pacific), -0.01, 0.01,max(Pacific)))) +
#   borders(database = "world2") + new_scale_color() + labs(x="", y = "") + 
#   coord_fixed(xlim = xlims%%360,ylim=ylims)  +geom_point(size = 1, col = "black", alpha = 0.3) + #scale_colour_distiller(palette = "Set1")+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(), legend.position = "none")# + facet_wrap(species_id ~ .)
# print(P)
# 
# 
# tiff('map_k0_chlocx.tiff', units="cm", width=10, height=16, res=800, compression = 'lzw')
# P
# dev.off()
# 
# # stop cluster use when we are done
# # parallel::stopCluster(cl = my.cluster)
# 
# save.image(file = "Synthetic_chl_ocx.RData")
# sendmail('jpinti@udel.edu', 'jpinti@udel.edu', 'Biased walks for chl ocx', '100 tracks wouhou!',control=list(smtpServer="mail.udel.edu"))


#################################
# Cleaning after everything has run
rm(maxsst, nc, newp, oceans, P, Simulation, Simulationstart, SST, a, B, comptbreak, concentration, durtrack, i, lat, latmini, lonmini,
   lon, maxpos, maxsstlat, maxsstlon, newd, pos, runnumber, satellitetime, StartD, Startlat, Startlon, STEP, Time, whichlat, whichlon, xlims, ylims, ANGLE, point.in.SpatialPolygons)

for (i in 1:100){
  TOT[[i]]$EV_ID <- paste(TOT[[i]]$kappa,TOT[[i]]$event_id,sep="_")
  TOT[[i]]$event_id <- as.numeric(TOT[[i]]$event_id)
  TOT[[i]]$kappa <- as.numeric(TOT[[i]]$kappa)
  TOT[[i]]$event_id <- as.numeric(TOT[[i]]$event_id)
  TOT[[i]]$longitude <- as.numeric(TOT[[i]]$longitude)
  TOT[[i]]$latitude <- as.numeric(TOT[[i]]$latitude)
  TOT[[i]]$bearing <- as.numeric(TOT[[i]]$bearing)
  TOT[[i]]$steplength <- as.numeric(TOT[[i]]$steplength)
  TOT[[i]]$dt <- as.numeric(TOT[[i]]$dt)
  
}

Synthetic <- TOT[[1]]
for (i in 2:100){
  Synthetic <- bind_rows(Synthetic, TOT[[i]])
  
}
rm(i,TOT)
save.image(file = "Fakeanimals_chl_ocx_0.RData")
