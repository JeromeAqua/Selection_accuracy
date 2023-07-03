# Merge all sub datasets together
library("glue")
library("dplyr")

VAR <- "ftle"
# setwd("/home/jpinti/Importance_accuracy/Extract_satellite_data")
# setwd("..")
load("../all_tracks0C.RData")
# setwd("For_0")

#For synthetic 1 
load("Accuracy_lcs_gaps_matched_2_kappa0_1.RData")
outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outT$sigma <- 1

outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})

# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB1 <- outB[[1]] 
outL1 <- outL[[1]] 
outC1 <- outC[[1]]
outJ1 <- outJ[[1]]
# outH1 <- outH[[1]]

for (i in 2:100){
  outB1 <- rbind(outB1, outB[[i]])
  outC1 <- rbind(outC1, outC[[i]])
  outL1 <- rbind(outL1, outL[[i]])
  outJ1 <- rbind(outJ1, outJ[[i]])
  # outH1 <- rbind(outH1, outH[[i]])
}


outB1$sigma <- 1
outC1$sigma <- 1
outL1$sigma <- 1
outJ1$sigma <- 1
# outH1$sigma <- 1

outT1 <- outT
outT1$event_id <- Synthetic$event_id
outT1$kappa <- Synthetic$kappa
outT1$date <- Synthetic$date
outT1$latitude <- Synthetic$latitude
outT1$longitude <- Synthetic$longitude
outT1$cue <- Synthetic$cue
outT1$EV_ID <- Synthetic$EV_ID


#For synthetic 10
load("Accuracy_lcs_gaps_matched_2_kappa0_10.RData")
outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outT$sigma <- 10

outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})

# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB10 <- outB[[1]] 
outL10 <- outL[[1]] 
outC10 <- outC[[1]]
outJ10 <- outJ[[1]]
# outH10 <- outH[[1]]

for (i in 2:100){
  outB10 <- rbind(outB10, outB[[i]])
  outC10 <- rbind(outC10, outC[[i]])
  outL10 <- rbind(outL10, outL[[i]])
  outJ10 <- rbind(outJ10, outJ[[i]])
  # outH10 <- rbind(outH10, outH[[i]])
}


outB10$sigma <- 10
outC10$sigma <- 10
outL10$sigma <- 10
outJ10$sigma <- 10
# outH10$sigma <- 10

outT10 <- outT
outT10$event_id <- Synthetic$event_id
outT10$kappa <- Synthetic$kappa
outT10$date <- Synthetic$date
outT10$latitude <- Synthetic$latitude
outT10$longitude <- Synthetic$longitude
outT10$cue <- Synthetic$cue
outT10$EV_ID <- Synthetic$EV_ID

#For synthetic 25
load("Accuracy_lcs_gaps_matched_2_kappa0_25.RData")
outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outT$sigma <- 25

outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})

# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB25 <- outB[[1]] 
outL25 <- outL[[1]] 
outC25 <- outC[[1]]
outJ25 <- outJ[[1]]
# outH25 <- outH[[1]]

for (i in 2:100){
  outB25 <- rbind(outB25, outB[[i]])
  outC25 <- rbind(outC25, outC[[i]])
  outL25 <- rbind(outL25, outL[[i]])
  outJ25 <- rbind(outJ25, outJ[[i]])
  # outH25 <- rbind(outH25, outH[[i]])
}


outB25$sigma <- 25
outC25$sigma <- 25
outL25$sigma <- 25
outJ25$sigma <- 25
# outH25$sigma <- 25

outT25 <- outT
outT25$event_id <- Synthetic$event_id
outT25$kappa <- Synthetic$kappa
outT25$date <- Synthetic$date
outT25$latitude <- Synthetic$latitude
outT25$longitude <- Synthetic$longitude
outT25$cue <- Synthetic$cue
outT25$EV_ID <- Synthetic$EV_ID


#For synthetic 50
load("Accuracy_lcs_gaps_matched_2_kappa0_50.RData")
outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outT$sigma <- 50

outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})

# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB50 <- outB[[1]] 
outL50 <- outL[[1]] 
outC50 <- outC[[1]]
outJ50 <- outJ[[1]]
# outH50 <- outH[[1]]

for (i in 2:100){
  outB50 <- rbind(outB50, outB[[i]])
  outC50 <- rbind(outC50, outC[[i]])
  outL50 <- rbind(outL50, outL[[i]])
  outJ50 <- rbind(outJ50, outJ[[i]])
  # outH50 <- rbind(outH50, outH[[i]])
}


outB50$sigma <- 50
outC50$sigma <- 50
outL50$sigma <- 50
outJ50$sigma <- 50
# outH50$sigma <- 50

outT50 <- outT
outT50$event_id <- Synthetic$event_id
outT50$kappa <- Synthetic$kappa
outT50$date <- Synthetic$date
outT50$latitude <- Synthetic$latitude
outT50$longitude <- Synthetic$longitude
outT50$cue <- Synthetic$cue
outT50$EV_ID <- Synthetic$EV_ID


#For synthetic 111
load("Accuracy_lcs_gaps_matched_2_kappa0_111.RData")
outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outT$sigma <- 111
outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})

# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB111 <- outB[[1]] 
outL111 <- outL[[1]] 
outC111 <- outC[[1]]
outJ111 <- outJ[[1]]
# outH111 <- outH[[1]]

for (i in 2:100){
  outB111 <- rbind(outB111, outB[[i]])
  outC111 <- rbind(outC111, outC[[i]])
  outL111 <- rbind(outL111, outL[[i]])
  outJ111 <- rbind(outJ111, outJ[[i]])
  # outH111 <- rbind(outH111, outH[[i]])
}

outB111$sigma <- 111
outC111$sigma <- 111
outL111$sigma <- 111
outJ111$sigma <- 111
# outH111$sigma <- 111

outT111 <- outT
outT111$event_id <- Synthetic$event_id
outT111$kappa <- Synthetic$kappa
outT111$date <- Synthetic$date
outT111$latitude <- Synthetic$latitude
outT111$longitude <- Synthetic$longitude
outT111$cue <- Synthetic$cue
outT111$EV_ID <- Synthetic$EV_ID

#For synthetic 0  #error here, correct name of files if we go further
load("Accuracy_lcs_gaps_matched_2_kappa0_1.RData")
outT$sigma <- 0

outT$WM_ftle <- outT$ftle #all this because 0 standard deviation
for (i in 1:100){
  outB[[i]]$WM_ftle <- outB[[i]]$ftle
  outL[[i]]$WM_ftle <- outL[[i]]$ftle
  outC[[i]]$WM_ftle <- outC[[i]]$ftle
  outJ[[i]]$WM_ftle <- outJ[[i]]$ftle
  # outH[[i]]$WM_sst <- outH[[i]]$sst
}

outT <- subset(outT, select = -c(ftle))
outB <- lapply(outB, function(DF){subset(DF, select = -c(ftle))})
outC <- lapply(outC, function(DF){subset(DF, select = -c(ftle))})
outL <- lapply(outL, function(DF){subset(DF, select = -c(ftle))})
outJ <- lapply(outJ, function(DF){subset(DF, select = -c(ftle))})

outB<- lapply(outB, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outB<- lapply(outB, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outB<- lapply(outB, function(DF){cbind(DF, date = Synthetic$date)})
outB<- lapply(outB, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outB<- lapply(outB, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outB<- lapply(outB, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outB<- lapply(outB, function(DF){cbind(DF, cue = Synthetic$cue)})

outL<- lapply(outL, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outL<- lapply(outL, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outL<- lapply(outL, function(DF){cbind(DF, date = Synthetic$date)})
outL<- lapply(outL, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outL<- lapply(outL, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outL<- lapply(outL, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outL<- lapply(outL, function(DF){cbind(DF, cue = Synthetic$cue)})

outC<- lapply(outC, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outC<- lapply(outC, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outC<- lapply(outC, function(DF){cbind(DF, date = Synthetic$date)})
outC<- lapply(outC, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outC<- lapply(outC, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outC<- lapply(outC, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outC<- lapply(outC, function(DF){cbind(DF, cue = Synthetic$cue)})

outJ<- lapply(outJ, function(DF){cbind(DF, event_id = Synthetic$event_id)})
outJ<- lapply(outJ, function(DF){cbind(DF, kappa = Synthetic$kappa)})
outJ<- lapply(outJ, function(DF){cbind(DF, date = Synthetic$date)})
outJ<- lapply(outJ, function(DF){cbind(DF, latitude = Synthetic$latitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, longitude = Synthetic$longitude)})
outJ<- lapply(outJ, function(DF){cbind(DF, EV_ID = Synthetic$EV_ID)})
outJ<- lapply(outJ, function(DF){cbind(DF, cue = Synthetic$cue)})
# outH<- lapply(outH, function(DF){cbind(DF, event_id = Synthetic$event_id)})
# outH<- lapply(outH, function(DF){cbind(DF, kappa = Synthetic$kappa)})
# outH<- lapply(outH, function(DF){cbind(DF, date = Synthetic$date)})
# outH<- lapply(outH, function(DF){cbind(DF, latitude = Synthetic$latitude)})
# outH<- lapply(outH, function(DF){cbind(DF, longitude = Synthetic$longitude)})

outB0 <- outB[[1]] 
outL0 <- outL[[1]] 
outC0 <- outC[[1]]
outJ0 <- outJ[[1]]
# outH0 <- outH[[1]]

for (i in 2:100){
  outB0 <- rbind(outB0, outB[[i]])
  outC0 <- rbind(outC0, outC[[i]])
  outL0 <- rbind(outL0, outL[[i]])
  outJ0 <- rbind(outJ0, outJ[[i]])
  # outH0 <- rbind(outH0, outH[[i]])
}


outB0$sigma <- 0
outC0$sigma <- 0
outL0$sigma <- 0
outJ0$sigma <- 0
# outH0$sigma <- 0

outT0 <- outT
outT0$event_id <- Synthetic$event_id
outT0$kappa <- Synthetic$kappa
outT0$date <- Synthetic$date
outT0$latitude <- Synthetic$latitude
outT0$longitude <- Synthetic$longitude
outT0$cue <- Synthetic$cue
outT0$EV_ID <- Synthetic$EV_ID

#Then merge everything together
outT <- rbind(outT0, outT1, outT10, outT25, outT50, outT111)
outB <- rbind(outB0, outB1, outB10, outB25, outB50, outB111)
outL <- rbind(outL0, outL1, outL10, outL25, outL50, outL111)
outC <- rbind(outC0, outC1, outC10, outC25, outC50, outC111)
outJ <- rbind(outJ0, outJ1, outJ10, outJ25, outJ50, outJ111)
# outH <- rbind(outH0, outH1, outH10, outH25, outH50, outH111)

outT$cue <- "ftle"
outB$cue <- "ftle"
outC$cue <- "ftle"
outL$cue <- "ftle"
outJ$cue <- "ftle"
save(outT, outB, outL, outC, outJ, file = "Synthetic_all_LCS_gaps_0.RData")
