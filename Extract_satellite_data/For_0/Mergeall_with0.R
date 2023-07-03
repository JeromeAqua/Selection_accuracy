#Merge all and for_0 datasets together
library(tidyverse)
setwd("/home/jpinti/Importance_accuracy/Extract_satellite_data")

#For LCS
load("Synthetic_all_LCS_0C.RData")
outT0 <- outT
outB0 <- outB
outC0 <- outC
outJ0 <- outJ
outL0 <- outL

load("Synthetic_LCS_onlylcs.RData")

outT <- subset(outT, kappa !=0)
outB <- subset(outB, kappa !=0)
outC <- subset(outC, kappa !=0)
outL <- subset(outL, kappa !=0)
outJ <- subset(outJ, kappa !=0)

outT <- rbind(outT0, outT)
outB <- rbind(outB0, outB)
outL <- rbind(outL0, outL)
outC <- rbind(outC0, outC)
outJ <- rbind(outJ0, outJ)

rm(outT0, outB0, outL0, outC0, outJ0)

save.image("Synthetic_LCS_onlylcs_v2.RData")

rm(list = ls())

# For chl
#For LCS
load("Synthetic_all_chl_0C.RData")
outT0 <- outT
outB0 <- outB
outC0 <- outC
outJ0 <- outJ
outL0 <- outL

load("Synthetic_chl_onlychl.RData")

outT <- subset(outT, kappa !=0)
outB <- subset(outB, kappa !=0)
outC <- subset(outC, kappa !=0)
outL <- subset(outL, kappa !=0)
outJ <- subset(outJ, kappa !=0)

outT <- rbind(outT0, outT)
outB <- rbind(outB0, outB)
outL <- rbind(outL0, outL)
outC <- rbind(outC0, outC)
outJ <- rbind(outJ0, outJ)

rm(outT0, outB0, outL0, outC0, outJ0)

save.image("Synthetic_chl_onlychl_v2.RData")

rm(list = ls())

# For SST
load("Synthetic_all_sst_0C.RData")
outT0 <- outT
outB0 <- outB
outC0 <- outC
outJ0 <- outJ
outL0 <- outL

load("Synthetic_SST_onlysst.RData")

outT <- subset(outT, kappa !=0)
outB <- subset(outB, kappa !=0)
outC <- subset(outC, kappa !=0)
outL <- subset(outL, kappa !=0)
outJ <- subset(outJ, kappa !=0)

outT <- rbind(outT0, outT)
outB <- rbind(outB0, outB)
outL <- rbind(outL0, outL)
outC <- rbind(outC0, outC)
outJ <- rbind(outJ0, outJ)

rm(outT0, outB0, outL0, outC0, outJ0)

save.image("Synthetic_SST_onlysst_v2.RData")