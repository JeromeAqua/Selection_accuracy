library("glue")
library("dplyr")

## Create dataframe for final matrices
#Fixing kappa, accuracy is y axis, # of data points is x axis 

#Here done for chlorophyll
k <- 2
CUE <- "chlorophyll"

Y <-  c(0,1,10,25,50,111) #Different accuracis

load(glue("Boostrapestimates_{CUE}_K{k}_v4.RData"))

ntracks <- unique(Allbootstrapestimates$ntracks)

# if (sigma == 0){ npointspertrack <- 32
# } else if (sigma == 1){ npointspertrack <- 48
# } else if (sigma == 10){ npointspertrack <- 69
# } else if (sigma == 25){ npointspertrack <- 77
# } else if (sigma %in% c(50,111)) {npointspertrack <- 80}


XAXIS <- 32:8000

Mless <- matrix(NA,nrow=6,ncol=length(XAXIS))
Mgreat <- Mless

for (s in 1:6){
  S <- Y[s]
  
  #First defining the resolution of data we have
  if (S == 0){ npointspertrack <- 32
  } else if (S == 1){ npointspertrack <- 48
  } else if (S == 10){ npointspertrack <- 69
  } else if (S == 25){ npointspertrack <- 77
  } else if (S %in% c(50,111)) {npointspertrack <- 80}
  
  for (x in 1:15){
    idx <- which(XAXIS==(ntracks[x]*npointspertrack))
    Mless[s,idx] <- Allbootstrapestimates$muAllless[Allbootstrapestimates$sigma==S&Allbootstrapestimates$ntracks==ntracks[x]]
    Mgreat[s,idx] <- Allbootstrapestimates$muAllgreat[Allbootstrapestimates$sigma==S&Allbootstrapestimates$ntracks==ntracks[x]]
  }
}

write.table(Mless, file = glue("Mless_{CUE}_{k}_v4.txt"), row.names = FALSE)
write.table(Mgreat, file = glue("Mgreat_{CUE}_{k}_v4.txt"), row.names = FALSE)
