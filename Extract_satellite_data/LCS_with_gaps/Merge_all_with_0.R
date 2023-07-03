# Merge all with 0

#Merge all and for_0 datasets together
library(tidyverse)
setwd("/home/jpinti/Importance_accuracy/Extract_satellite_data/LCS_with_gaps")

#For LCS
load("Synthetic_all_LCS_gaps_0.RData")
outT0 <- outT
outB0 <- outB
outC0 <- outC
outJ0 <- outJ
outL0 <- outL

load("Synthetic_all_LCS_gaps.RData")

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

save.image("Synthetic_LCSgaps_onlylcs.RData")