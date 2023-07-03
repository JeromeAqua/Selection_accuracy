#Merge the three datasets
load("Fakeanimals_chl_ocx.RData")
Syntheticchl <- Synthetic
Syntheticchl$cue <- "chlorophyll"
Syntheticchl$EV_ID <- paste(Syntheticchl$EV_ID,"chl",sep="_")

load("Fakeanimals_SST.RData")
Syntheticsst <- Synthetic
Syntheticsst$cue <- "sst"
Syntheticsst$EV_ID <- paste(Syntheticsst$EV_ID,"sst",sep="_")

load("Fakeanimals_ftle.RData")
Syntheticftle <- Synthetic
Syntheticftle$cue <- "ftle"
Syntheticftle$EV_ID <- paste(Syntheticftle$EV_ID,"ftle",sep="_")


Synthetic <- rbind(Syntheticsst, Syntheticchl, Syntheticftle)

save(Synthetic, file = "Fakeanimals.RData")
