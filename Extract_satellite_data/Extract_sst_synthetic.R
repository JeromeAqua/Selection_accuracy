args <- commandArgs(TRUE)
# J <- strtoi(args[1], base=10L)


jj <- strtoi(args[1], base=10L)#floor(J[1]/4) + 1 # This is the variable we work with -- 37 for sst
# 
# tt <- (J[1] %% 4) 

iii <- 0 #standard deviation for the uncertainty

print(iii)

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
load("all_tracks0C.RData")
Trackingdates <- unique(Synthetic$date)
Trackingdates <- ymd_hms(Trackingdates)
# Trackingdates <- Trackingdates[as.numeric(format(Trackingdates, format = "%Y")) >2001] ## Remove data points before 2001 as we won't be able to interpolate anyways

# Limits$XMIN[LIMITS$XMIN < -180] <- -180
# Limits$XMAX[LIMITS$XMAX > 180] <- 180

## Load satellite data
ModisAquaDaily <- nc_open("/lustre/scratch/orb/aqua/8day/aqua_8day.2006.concat.nc4", write = FALSE)
yr <- 2006
#"http://basin.ceoe.udel.edu/thredds/dodsC/AQUAGLOBAL8DAY9KM.nc")
# setwd("/home/2453/Extract_satellite_data/")

Latitude <- ncvar_get(ModisAquaDaily, varid = "lat", start = 1, count = -1)
Longitude <- ncvar_get(ModisAquaDaily, varid = "lon", start = 1, count = -1)
Time <- ncvar_get(ModisAquaDaily, varid = "time", start = 1, count = -1)
Time <- as.POSIXct(Time, origin = "1970-01-01", tz = "UTC")

Products <- names(ModisAquaDaily$var)
nlon <- length(Longitude)

Products
print(glue("{Products[jj]}"))

## Initialize the datasets with NA columns for the ocean products we want  ("+1" is here to have class = numerical and not logical)

Synthetic <- ungroup(Synthetic)  #somehow needed to make cbind work on caviness
Synthetic <- cbind(Synthetic, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj]) )
Synthetic <- cbind(Synthetic, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}")) )
Synthetic <- cbind(Synthetic, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}")) )

BR <- lapply(BR, function(DF) {cbind(DF, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj])) })
BR <- lapply(BR, function(DF) {cbind(DF, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}"))) })
BR <- lapply(BR, function(DF) {cbind(DF, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}"))) })
BR <- lapply(BR, function(DF){cbind(DF, setNames( lapply(glue("id"), function(x) x=seq.int(nrow(DF))), glue("id"))) } )

LF <- lapply(LF, function(DF) {cbind(DF, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj])) })
LF <- lapply(LF, function(DF) {cbind(DF, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}"))) })
LF <- lapply(LF, function(DF) {cbind(DF, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}"))) })
LF <- lapply(LF, function(DF){cbind(DF, setNames( lapply(glue("id"), function(x) x=seq.int(nrow(DF))), glue("id"))) } )

CRW <- lapply(CRW, function(DF) {cbind(DF, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj])) })
CRW <- lapply(CRW, function(DF) {cbind(DF, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}"))) })
CRW <- lapply(CRW, function(DF) {cbind(DF, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}"))) })
CRW <- lapply(CRW, function(DF){cbind(DF, setNames( lapply(glue("id"), function(x) x=seq.int(nrow(DF))), glue("id"))) } )

JD <- lapply(JD, function(DF) {cbind(DF, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj])) })
JD <- lapply(JD, function(DF) {cbind(DF, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}"))) })
JD <- lapply(JD, function(DF) {cbind(DF, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}"))) })
JD <- lapply(JD, function(DF){cbind(DF, setNames( lapply(glue("id"), function(x) x=seq.int(nrow(DF))), glue("id"))) } )

# 
# for (hhh in 1:100){
#   Hull[[hhh]] <- ungroup(Hull[[hhh]])
#   Hull$longitude[Hull$longitude > 180] <- Hull$longitude[Hull$longitude > 180]-360
# }
# Hull <- lapply(Hull, function(DF) {cbind(DF, setNames( lapply(Products[jj], function(x) x=NA+1), Products[jj])) })
# Hull <- lapply(Hull, function(DF) {cbind(DF, setNames( lapply(glue("WM_{Products[jj]}"), function(x) x=NA+1), glue("WM_{Products[jj]}"))) })
# Hull <- lapply(Hull, function(DF) {cbind(DF, setNames( lapply(glue("SD_{Products[jj]}"), function(x) x=NA+1), glue("SD_{Products[jj]}"))) })

#Put the IQR in all datasets -- easier to deal with / sigma = IQR/4
IQR <- iii*4* 0.008 #in degrees

Synthetic$LonErrorIQR <- IQR
Synthetic$LatErrorIQR <- IQR
for (dst in 1:100){
  BR[[dst]]$LonErrorIQR <- IQR
  BR[[dst]]$LatErrorIQR <- IQR

  LF[[dst]]$LonErrorIQR <- IQR
  LF[[dst]]$LatErrorIQR <- IQR
  
  CRW[[dst]]$LonErrorIQR <- IQR
  CRW[[dst]]$LatErrorIQR <- IQR
  
  JD[[dst]]$LonErrorIQR <- IQR
  JD[[dst]]$LatErrorIQR <- IQR
}

## Function to generate normally distributed values based on N((mux, muy), (sx,sy))
gaussian = function(x,y,mux,muy,sx,sy){
  return(exp(-1/2 * ( ((x %% 360 -mux %% 360)/sx)^2 + ((y-muy)/sy)^2 )  ))
}

## Function to compute the weighted mean and standard deviation of Field around (mux,muy) with interquartile range (ix,iy)
WM_and_SD = function(mux,muy,ix,iy,Field){
  if (!is.na(mux) & !is.na(muy)){
    #Make longitude vector - careful to wrap correctly when crossing the date line
    Minlon <- which.min(abs((Longitude %% 360) - ((mux-ix) %%360) ))
    Maxlon <- which.min(abs((Longitude %% 360) - ((mux+ix) %%360) ))
    if (Minlon <= Maxlon){
      longs <- Longitude[Minlon:Maxlon]
      lonidx <- Minlon:Maxlon
    } else{
      longs <- c(Longitude[1:Maxlon], Longitude[Minlon:Maxlon])
      lonidx <- c(1:Maxlon, Minlon:nlon)
    }
    
    #Make latitude vector 
    Minlat <- which.min(abs(Latitude - (muy - iy)))
    Maxlat <- which.min(abs(Latitude - (muy + iy)))
    lats <- Latitude[Minlat:Maxlat]
    
    #Make meshgrid for interpolation
    MG <- meshgrid(longs, lats)
    MGlongs <- c(MG$X)
    MGlats <- c(MG$Y)
    
    Matrixfield <- Field[lonidx, Minlat:Maxlat]
    variablevalues <- c(Matrixfield)
    
    #STD is the standard deviation of the values in the confidence interval
    STD <- sd(variablevalues, na.rm = TRUE)
    
    # Now compute the weighted mean of the area
    Test <- !is.na(variablevalues)
    
    if (sum(Test) > 0 ){
    #  ############################################################################################################################################################ MODIFIED HERE
      Matgaussfactors <- mapply(function(x,y) {return(gaussian(x,y,mux,muy,ix/4,iy/4))}, MGlongs, MGlats)
      if (sum(Matgaussfactors)==0) { Matgaussfactors <- 1}
      Wmean <- sum(Matgaussfactors*variablevalues / sum(Matgaussfactors*Test), na.rm = TRUE)
      # Wmean <- mean(variablevalues, na.rm = TRUE)
    } else {
      Wmean <- NA
    }
    
    return(c(Wmean, STD))
  }
  else {
    return(c(NA,NA))
  }
}

## Function to interpolate easily the values of nameprod from Field at corresponding values in DF       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!No WM and SD for synthetic tracks
Extract = function(DF, Field, nameprod, ids){
  Values <- interp2(Latitude, Longitude, Field, DF$latitude[DF$id %in% ids], DF$longitude[DF$id %in% ids])
  
  DF[[nameprod]][DF$id %in% ids] <- Values #DF$id %in% 
  
  W_S <- c()
  for (ii in ids){
    W_S <- rbind(W_S, WM_and_SD(DF$longitude[ii], DF$latitude[ii], DF$LonErrorIQR[ii], DF$LatErrorIQR[ii], Field))
  }

  DF[[glue("WM_{Products[jj]}")]][DF$id %in% ids] <- W_S[,1]
  DF[[glue("SD_{Products[jj]}")]][DF$id %in% ids] <- W_S[,2]
  
  return(DF)
}


##################################################################################
######################### EXTRACTION STARTS HERE #################################
##################################################################################


for(i in 1:length(Trackingdates)){ #Start at Trackingdates[103] because no satellite data before 2002
  # as.numeric(format(Trackingdates[1000], format = "%Y"))
  print(i)
  # strftime(Trackingdates[1000], format = "%j")
  
  tryCatch({
    ## Here close and open new one if we change year
    if ( as.numeric(format(Trackingdates[i], format = "%Y")) != yr  ){
      nc_close(ModisAquaDaily)
      yr <- as.numeric(format(Trackingdates[i], format = "%Y"))
      filetoopen <- as.character(glue("/lustre/scratch/orb/aqua/8day/aqua_8day.{yr}.concat.nc4"))
      ModisAquaDaily <- nc_open(filetoopen, write = FALSE)
      Time <- ncvar_get(ModisAquaDaily, varid = "time", start = 1, count = -1)
      Time <- as.POSIXct(Time, origin = "1970-01-01", tz = "UTC")
    }
    
    
    if (min(abs(difftime(Time, Trackingdates[i], unit = "days"))) < 1){ #only interpolate if there is actual data from that day
      
      satellitetime <- which.min(abs(Time - Trackingdates[i])) #the satellite time we are going to interpolate from
      
      ids <- Synthetic$id[ymd_hms(Synthetic$date) == Trackingdates[i]]
      
      compt <- 0
      success <- 0
      while ( (success==0) & (compt<10) ){
        
        tryCatch(
          {
            Field <- ncvar_get(ModisAquaDaily, varid = Products[jj], start = c(1,1,satellitetime), count = c(-1,-1,1))  # then we extract what we need from the satellites
            
            
            # a <-proc.time()
            Synthetic <- Extract(Synthetic, Field, Products[jj], ids)
            
            
            # Brownian = foreach(b = 1:100, .export = c("Extract", "gaussian", "WM_and_SD"), .packages =
            #                     c("tidyverse", "glue", "dplyr","ncdf4","data.table", "pracma"), .combine = "c" ) %dopar%
            #  {
            #   outy <- Extract(Brownian[[b]], Field, Products[jj], ids)
            #   list(outy)
            # }
            BR <- mclapply(BR, function(DF) {Extract(DF, Field, Products[jj], ids)}, mc.cores = 10)
            LF <-     mclapply(LF,     function(DF) {Extract(DF, Field, Products[jj], ids)}, mc.cores = 10)
            CRW <-      mclapply(CRW,      function(DF) {Extract(DF, Field, Products[jj], ids)}, mc.cores = 10)
            JD <-       mclapply(JD,       function(DF) {Extract(DF, Field, Products[jj], ids)}, mc.cores = 10)
            # Hull <-     mclapply(Hull,     function(DF) {Extract(DF, Field, Products[jj], ids)}, mc.cores = 10)
            
            # print(proc.time() - a)
            
            success <- 1
          },
          error = function(e){
            print(glue("{i}: Didnt work!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
            print(compt)
          },
          finally = {compt <- compt + 1}
        )
      }
      
      if (compt == 10){
        # save(Brownian, Levy, CRW, JD, Synthetic, file = paste(glue("B_1day_9km_variable_{jj}.RData")))
        print("COULD NOT GET THIS DAY TO WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }
      
      if (i %%100 ==0 ){
        print(glue("{i} / 3116"))
      }
    }
    
  },
  error=function(e) {
    print("Could not extract data for this day")  
    return(i)
  } )
}

# save.image(paste(glue("/lustre/scratch/jpinti/Debug_{jj}.RData")))
## Now make all extracted data in a friendly format to facilitate saves and future work
outT <- data.frame(id = Synthetic$id)
# print("229")
outT[Products[jj]] = Synthetic[Products[jj]]
# print("230")
outT[glue("WM_{Products[jj]}")] <- Synthetic[glue("WM_{Products[jj]}")]
outT[glue("SD_{Products[jj]}")] <- Synthetic[glue("SD_{Products[jj]}")]
# print("235")

outB <- list()
for (a in 1:100){
  K <- data.frame(id = BR[[a]]$id)
  # if (a==1){print("240")}
  K[Products[jj]] = BR[[a]][Products[jj]]
  K[glue("WM_{Products[jj]}")] <- BR[[a]][glue("WM_{Products[jj]}")]
  K[glue("SD_{Products[jj]}")] <- BR[[a]][glue("SD_{Products[jj]}")]
  outB <- c(outB, list(K))
}

outL <- list()
for (a in 1:100){
  K <- data.frame(id = LF[[a]]$id)
  # if (a==1){print("250")}
  K[Products[jj]] = LF[[a]][Products[jj]]
  K[glue("WM_{Products[jj]}")] <- LF[[a]][glue("WM_{Products[jj]}")]
  K[glue("SD_{Products[jj]}")] <- LF[[a]][glue("SD_{Products[jj]}")]
  outL <- c(outL, list(K))
}

outC <- list()
for (a in 1:100){
  K <- data.frame(id = CRW[[a]]$id)
  K[Products[jj]] = CRW[[a]][Products[jj]]
  K[glue("WM_{Products[jj]}")] <- CRW[[a]][glue("WM_{Products[jj]}")]
  K[glue("SD_{Products[jj]}")] <- CRW[[a]][glue("SD_{Products[jj]}")]
  
  outC <- c(outC, list(K))
}

outJ <- list()
for (a in 1:100){
  K <- data.frame(id = JD[[a]]$id)
  K[Products[jj]] = JD[[a]][Products[jj]]
  K[glue("WM_{Products[jj]}")] <- JD[[a]][glue("WM_{Products[jj]}")]
  K[glue("SD_{Products[jj]}")] <- JD[[a]][glue("SD_{Products[jj]}")]
  
  outJ <- c(outJ, list(K))
}

# outH <- list()
# for (a in 1:100){
#   K <- data.frame(id = Hull[[a]]$id)
#   K[Products[jj]] = Hull[[a]][Products[jj]]
#   K[glue("WM_{Products[jj]}")] <- Hull[[a]][glue("WM_{Products[jj]}")]
#   K[glue("SD_{Products[jj]}")] <- Hull[[a]][glue("SD_{Products[jj]}")]
#   
#   outH <- c(outH, list(K))
# }


nc_close(ModisAquaDaily)
warnings()

# #stop cluster use when we are done
# parallel::stopCluster(cl = my.cluster)


save(outT, outB, outC, outL, outJ, file = paste(glue("/lustre/scratch/jpinti/Accuracy/Accuracy_matched_{jj}_{iii}_0C.RData")))

# save.image("Seawifs_try.RData")
# sendmail('jpinti@udel.edu', 'jpinti@udel.edu', 'Try of seawifs extraction over', 'Fingers crossed!',control=list(smtpServer="mail.udel.edu"))
# 
# #################################################################################
# #################################################################################
# ####################### EXTRACT VALUES FOR CONVEX HULL ##########################
# #################################################################################
# #################################################################################
# library("sp")
# library("raster")
# library("lubridate")
# library("rgeos")
# library("dplyr") #loaded after to remove a bug
# ###############################################################################
# # Function to know if a point is in the ocean
# point.in.SpatialPolygons = function(point.x, point.y, SpP){
#   # VERBATIM COPY FROM prevR package. They deserve all the credit for this function
#   ###############################################################################################
#   # Cette fonction renvoie pour chaque point defini par le couple (point.x, point.y) T ou F 
#   #     si le point est a l'interieur ou non du spatialPolygons SpP
#   # Un point est considere a l'interieur de N polygons si il est a l'interieur d'au moins
#   # un polygon non Hole et a l'exterieur de tous les polygons Hole
#   # Cette foncion est utilisee par toutes les fonctions de lissage (krige , kde, idw) . 
#   # En effet ces fonctions travaillent sur un grid rectangulaire englobant les donnees. En presentation on ne veut que les resulats 
#   #   interieurs a la frontiere qui est definie dans l'element SpP (SpP contient boundary). 
#   #   Tous les elements du grid hors de la frontiere seront dans les programmes de lissage positonnes a NA
#   # 
#   ###############################################################################################
#   
#   X = slot(SpP,"polygons")
#   is.inside = F
#   for(i in 1:length(X)){
#     PS   = slot(X[[i]],"Polygons")
#     for(j in 1:length(PS)){
#       pol.x = slot(PS[[j]],"coords")[,1]
#       pol.y = slot(PS[[j]],"coords")[,2]
#       pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
#       if(!slot(PS[[j]],"hole")) {
#         is.inside = is.inside | pointsPosition != 0
#       }
#     }
#   }
#   is.outsideHole = T
#   for(i in 1:length(X)){
#     PS   = slot(X[[i]],"Polygons")
#     for(j in 1:length(PS)){
#       pol.x = slot(PS[[j]],"coords")[,1]
#       pol.y = slot(PS[[j]],"coords")[,2]
#       pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
#       if(slot(PS[[j]],"hole")) {
#         is.outsideHole = is.outsideHole & (pointsPosition == 0 |  pointsPosition == 3)
#       }
#     }
#   }
#   is.inside & is.outsideHole
# }
# #################################################################################
# 
# 
# 
# #Generate data frame with all possible longitudes and latitudes
# library("pracma")
# POS <- meshgrid(1:length(Longitude), 1:length(Latitude))
# POS <- data.frame(longidx = c(POS$X), latidx = c(POS$Y))
# Trackingdates <- unique(Synthetic$date)
# load("../Reanalysis_convex.RData")
# 
# for (i in 1:length(Trackingdates)){
#   
#   #Open the file we need
#   satellitetime <- which.min(abs(Time - Trackingdates[i]))
#   Field <- ncvar_get(ModisAquaDaily, varid = Products[jj], start = c(1,1,satellitetime), count = c(-1,-1,1))  # then we extract what we need from the satellites
#   
#   speciesthatday <- unique(TOPP$species_id[TOPP$date==Trackingdates[i]])
#   for (j in 1:length(speciesthatday)){
#     # print(Trackingdates[i])
#     # print(j)
#     #Make the convex hull
#     A<- na.omit(TOPP[TOPP$species_id==speciesthatday[j] & month(TOPP$date)==month(Trackingdates[i]),c("latitude","longitude")])
#     if (TOPP$species_id[j]!="51b"){ #For Atlantic Bluefin tuna we cross the 0 line so need to have long from -180 to 180
#       #Create the spatial polygon
#       ch <- chull(A$longitude%%360, A$latitude)
#       coords <- A[c(ch, ch[1]),]
#       coords <- coords[, (names(coords) %in% c("longitude","latitude"))]
#       coords <- coords[, c(2, 1)]
#       coords$longitude <- coords$longitude%%360
#     }
#     else
#     {
#       ch <- chull(A$longitude, A$latitude)
#       coords <- A[c(ch, ch[1]),]
#       coords <- coords[, (names(coords) %in% c("longitude","latitude"))]
#       coords <- coords[, c(2, 1)]
#     }
#     
#     sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
#     # Expand the spatial polygon (= buffer) for x(=1 here) degree(s) around it
#     p <- raster::buffer(sp_poly, width = 1)
#     
#     #Find the points in the convex hull
#     if (TOPP$species_id[j]!="51b"){
#       inside <- point.in.SpatialPolygons(Longitude[POS$longidx]%%360, Latitude[POS$latidx], p)
#     } else {
#       inside <- point.in.SpatialPolygons(Longitude[POS$longidx], Latitude[POS$latidx], p)  
#     }
#     
#     # #Take 1000 random samples if too big
#     if (sum(inside)>1000){
#       idinside <- 1:length(inside)
#       idinside <- idinside[inside]
#       idinside <- sample(idinside, 1000)
#       inside <- logical(length(inside))
#       inside[idinside] <- TRUE
#     }
# 
#     lons <- POS$longidx[inside]
#     lats <- POS$latidx[inside]
#     
#     #Get the variable values in the convex hull
#     values <- mapply(function(x,y) {Field[x,y]}, lons, lats)
#     
#     #Collate the values in outH
#     temp <- data.frame(matrix(ncol = 5, nrow = 0))
#     colnames(temp) <- cols
#     
#     temp <- data.frame(species_id = speciesthatday[j], date = Trackingdates[i], latitude = Latitude[lats], longitude = Longitude[lons], dumdum = values)
#     names(temp)[names(temp) == 'dumdum'] <- glue("WM_{Products[jj]}")
#     
#     outH <- rbind(outH, temp)
#     
#   }
#   print(i)
# }
# 
# save(outH, file = paste(glue("/lustre/scratch/jpinti/9km_1day/Modis_1day_9km_variable_{jj}H.RData")))
