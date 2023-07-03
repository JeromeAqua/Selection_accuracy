args <- commandArgs(TRUE)
J <- strtoi(args[1], base=10L)

# Do KS tests
setwd("/home/2453/Importance_accuracy/KS_tests")
library("ggplot2")
library("gridExtra")
library("glue")
# library("patchwork")
library("boot")

load("/lustre/scratch/jpinti/Accuracy/Synthetic_chl_onlychl_v2.RData")
# load("/lustre/scratch/jpinti/Accuracy/Synthetic_all_chl_0.RData")
CUE <- "chlorophyll"


outT <- outT[outT$cue==CUE,]
outB <- outB[outB$cue==CUE,]
outL <- outL[outL$cue==CUE,]
outC <- outC[outC$cue==CUE,]
outJ <- outJ[outJ$cue==CUE,]

#Create dataframe for all results
Colnames <- c("kappa","sigma","ntracks","muBrownless", "muLevyless", "muCRWless", "muJDless", "muAllless", "muBrowngreat", "muLevygreat", "muCRWgreat", "muJDgreat", "muAllgreat",
              "sdBrownless", "sdLevyless", "sdCRWless", "sdJDless", "sdAllless", "sdBrowngreat", "sdLevygreat", "sdCRWgreat", "sdJDgreat", "sdAllgreat")
Allbootstrapestimates <- setNames(data.frame(matrix(ncol = length(Colnames), nrow = 0)), Colnames)

Colnames <- c("kappa", "sigma", "ntracks", "pvalBrownLess",  "pvalLevyLess",  "pvalCRWLess",  "pvalJDLess",  "pvalallLess",
              "pvalBrownGreat", "pvalLevyGreat", "pvalCRWGreat", "pvalJDGreat", "pvalallGreat",
              "DBrownLess",  "DLevyLess",  "DCRWLess",  "DJDLess",  "DallLess",
              "DBrownGreat", "DLevyGreat", "DCRWGreat", "DJDGreat", "DallGreat")
Alltesttracks <- setNames(data.frame(matrix(ncol = length(Colnames), nrow = 0)), Colnames)

kappavector <- c(0.25, 0, 0.75, 2, 0.5,1,5,10,20) # 10 5 2 1 0.75 0.5 0.25 0 ##HERE WE DO IT ONE BY ONE ON BIOHEN TO HAVE THE RESULTS FASTER
kappavector <- kappavector[J[1]]
sigmavector <- c(0,1,10,25,50,111)

for (K in kappavector){
  gc()
  print("kappa")
  print(K)
  for (S in sigmavector){
    print("sigma")
    print(S)
    
    # ## Fix kappa and sigma and check how many tracks we need to get significant results - to do only one run
    # K <- 0
    # S <- 0
    
    Colnames <- c("kappa", "sigma", "ntracks", "pvalBrownLess",  "pvalLevyLess",  "pvalCRWLess",  "pvalJDLess",  "pvalallLess",
                  "pvalBrownGreat", "pvalLevyGreat", "pvalCRWGreat", "pvalJDGreat", "pvalallGreat",
                  "DBrownLess",  "DLevyLess",  "DCRWLess",  "DJDLess",  "DallLess",
                  "DBrownGreat", "DLevyGreat", "DCRWGreat", "DJDGreat", "DallGreat")
    
    testtracks <- setNames(data.frame(matrix(ncol = length(Colnames), nrow = 0)), Colnames)
    
    N <-100 #number of repeat per (number of tracks taken together)
    NTRACKS <- c(1,2,3,5,10,15,20,30,40,50,60,70,80,90,99)#set of number of tracks investigated
    
    
    for (ntracks in NTRACKS){
      thresh <- 0.05 #threshold for significance (if /100 -> ######################### CHANGED HERE BECAUSE WE TEST A BONFERRONI CORRECTION)
      Dthresh <- 0.05
      compt <- 0
      repeat { #Repeat loop so that we can do it 50 times easily
        
        #First we select ntracks events 
        Eventswewant <- sample(1:100, ntracks)
        # Ids <- outT$event_id[outT$kappa == K & outT$sigma==sigma]
        
        #Select the corresponding WM_chl_ocx values
        TRACKS <-   subset(outT,outT$event_id %in% Eventswewant & outT$kappa == K & outT$sigma==S)
        Brownian <- subset(outB,outB$event_id %in% Eventswewant & outB$kappa == K & outB$sigma==S) 
        Levy <-     subset(outL,outL$event_id %in% Eventswewant & outL$kappa == K & outL$sigma==S) 
        CRW <-      subset(outC,outC$event_id %in% Eventswewant & outC$kappa == K & outC$sigma==S) 
        JD <-       subset(outJ,outJ$event_id %in% Eventswewant & outJ$kappa == K & outJ$sigma==S) 
        
        # #Put all the tracks in a single dataframe -- no need they already are
        # Brownian <- bind_rows(Brownian, .id = "column_label") #column_label is the run# for all the synthetic tracks
        # Levy <- bind_rows(Levy, .id = "column_label")
        # CRW <- bind_rows(CRW, .id = "column_label")
        # JD <- bind_rows(JD, .id = "column_label")
        
        #Now do the KS tests
        ALT <- "less" # c("less","greater","two.sided")
        k <- ks.test(TRACKS$WM_chl_ocx,Brownian$WM_chl_ocx, alternative = ALT)
        D1 <- k$statistic
        p1 <- k$p.value
        
        # k <- ks.test(TRACKS$WM_chl_ocx,Levy$WM_chl_ocx, alternative = ALT)
        D2 <-  k$statistic
        p2 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,CRW$WM_chl_ocx, alternative = ALT)
        D3 <- k$statistic
        p3 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,JD$WM_chl_ocx, alternative = ALT)
        D4 <- k$statistic
        p4 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,c(JD$WM_chl_ocx,CRW$WM_chl_ocx,Brownian$WM_chl_ocx), alternative = ALT)
        D5 <- k$statistic
        p5 <- k$p.value
        
        ALT <- "greater" # c("less","greater","two.sided")
        k <- ks.test(TRACKS$WM_chl_ocx,Brownian$WM_chl_ocx, alternative = ALT)
        D6 <- k$statistic
        p6 <- k$p.value
        
        # k <- ks.test(TRACKS$WM_chl_ocx,Levy$WM_chl_ocx, alternative = ALT)
        D7 <- k$statistic
        p7 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,CRW$WM_chl_ocx, alternative = ALT)
        D8 <- k$statistic
        p8 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,JD$WM_chl_ocx, alternative = ALT)
        D9 <- k$statistic
        p9 <- k$p.value
        
        k <- ks.test(TRACKS$WM_chl_ocx,c(JD$WM_chl_ocx,CRW$WM_chl_ocxBrownian$WM_chl_ocx), alternative = ALT)
        D10 <- k$statistic
        p10 <- k$p.value
        
        
        testtracks[nrow(testtracks)+1,] <- c(K, S, ntracks, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                                             D1, D2, D3, D4, D5, D6, D7, D8, D9, D10)
        
        compt <- compt + 1
        if (compt == N){
          break
        }
      }
      print(ntracks)
    }
    
    
    #####################################################################################################
    # BOOTSTRAP ESTIMATE OF MEAN NUMBER OF SIGNIFICANT TEST PER NUMBER OF TRACKS LUMPED TOGETHER
    
    Colnames <- c("kappa","sigma","ntracks","muBrownless", "muLevyless", "muCRWless", "muJDless", "muAllless", "muBrowngreat", "muLevygreat", "muCRWgreat", "muJDgreat", "muAllgreat",
                  "sdBrownless", "sdLevyless", "sdCRWless", "sdJDless", "sdAllless", "sdBrowngreat", "sdLevygreat", "sdCRWgreat", "sdJDgreat", "sdAllgreat")
    ntrackBS <- setNames(data.frame(matrix(ncol = length(Colnames), nrow = length(NTRACKS))), Colnames)
    
    ntrackBS$kappa <- K
    ntrackBS$sigma <- S
    
    ntrackBS$ntracks <- NTRACKS
    ntrackBS[is.na(ntrackBS)==TRUE] <- 0
    
    for (nn in 1:length(NTRACKS)){
      ntr <- NTRACKS[nn]
      dftemp <- testtracks[testtracks$ntracks==ntr,]
      
      f1 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalBrownLess < thresh & d2$DBrownLess>Dthresh)/dim(d2)[1])
      }
      
      BS1 <- boot(dftemp, f1, R = 10000)
      ntrackBS$sdBrownless[nn] <- sd(BS1$t)
      ntrackBS$muBrownless[nn] <- mean(BS1$t)
      
      f2 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalBrownGreat < thresh & d2$DBrownGreat>Dthresh)/dim(d2)[1])
      }
      
      BS2 <- boot(dftemp, f2, R = 10000)
      ntrackBS$sdBrowngreat[nn] <- sd(BS2$t)
      ntrackBS$muBrowngreat[nn] <- mean(BS2$t)
      
      
      # f3 <- function(DF, i){
      #   d2 <- DF[i,]
      #   return(sum(d2$pvalLevyGreat < thresh & d2$DLevyGreat>Dthresh)/dim(d2)[1])
      # }
      
      # BS3 <- boot(dftemp, f3, R = 10000)
      # ntrackBS$sdLevygreat[nn] <- sd(BS3$t)
      # ntrackBS$muLevygreat[nn] <- mean(BS3$t)
      
      # f4 <- function(DF, i){
      #   d2 <- DF[i,]
      #   return(sum(d2$pvalLevyLess < thresh & d2$DLevyLess>Dthresh)/dim(d2)[1])
      # }
      
      # BS4 <- boot(dftemp, f4, R = 10000)
      # ntrackBS$sdLevyless[nn] <- sd(BS4$t)
      # ntrackBS$muLevyless[nn] <- mean(BS4$t)
      
      f5 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalCRWLess < thresh & d2$DCRWLess>Dthresh)/dim(d2)[1])
      }
      
      BS5 <- boot(dftemp, f5, R = 10000)
      ntrackBS$sdCRWless[nn] <- sd(BS5$t)
      ntrackBS$muCRWless[nn] <- mean(BS5$t)
      
      f6 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalCRWGreat < thresh & d2$DCRWGreat>Dthresh)/dim(d2)[1])
      }
      
      BS6 <- boot(dftemp, f6, R = 10000)
      ntrackBS$sdCRWgreat[nn] <- sd(BS6$t)
      ntrackBS$muCRWgreat[nn] <- mean(BS6$t)
      
      f7 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalJDGreat < thresh & d2$DJDGreat>Dthresh)/dim(d2)[1])
      }
      
      BS7 <- boot(dftemp, f7, R = 10000)
      ntrackBS$sdJDgreat[nn] <- sd(BS7$t)
      ntrackBS$muJDgreat[nn] <- mean(BS7$t)
      
      f8 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalJDLess < thresh & d2$DJDLess>Dthresh)/dim(d2)[1])
      }
      
      BS8 <- boot(dftemp, f8, R = 10000)
      ntrackBS$sdJDless[nn] <- sd(BS8$t)
      ntrackBS$muJDless[nn] <- mean(BS8$t)
      
      f9 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalallLess < thresh & d2$DallLess>Dthresh)/dim(d2)[1])
      }
      
      BS9 <- boot(dftemp, f9, R = 10000)
      ntrackBS$sdAllless[nn] <- sd(BS9$t)
      ntrackBS$muAllless[nn] <- mean(BS9$t)
      
      f10 <- function(DF, i){
        d2 <- DF[i,]
        return(sum(d2$pvalallGreat < thresh & d2$DallGreat>Dthresh)/dim(d2)[1])
      }
      
      BS10 <- boot(dftemp, f10, R = 10000)
      ntrackBS$sdAllgreat[nn] <- sd(BS10$t)
      ntrackBS$muAllgreat[nn] <- mean(BS10$t)
      
      print(nn)
    }
    ####################################################################################################
    
    # colors = c("black","red","blue","darkgoldenrod2","purple","chartreuse3")
    # szl <- 1
    # szp <- 2
    
   # Qabove <- ggplot(data = ntrackBS, aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), labels = c("","","",""),limits = c(1, NA)) +
   #   scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6,0.8,1)) +
   #   geom_ribbon(aes(ymin = pmax(0,muAllless-3*sdAllless), ymax = pmin(1,muAllless+3*sdAllless)), fill = colors[1], alpha = 0.2) +
   #   geom_ribbon(aes(ymin = pmax(0,muBrownless-3*sdBrownless), ymax = pmin(1,muBrownless + 3*sdBrownless)), fill = colors[2], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muLevyless-3*sdLevyless), ymax = pmin(1,muLevyless+3*sdLevyless)), fill = colors[3], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muCRWless-3*sdCRWless), ymax = pmin(1,muCRWless+3*sdCRWless)), fill = colors[5], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muJDless-3*sdJDless), ymax = pmin(1,muJDless+3*sdJDless)), fill = colors[6], alpha = 0.2) +
   #  # geom_ribbon(aes(ymin = Allgreatmed-Allgreatsd, ymax = Allgreatmed+Allgreatsd), fill = colors[1], alpha = 0.2) +
   #  # geom_ribbon(aes(ymin = Browngreatmed-Browngreatsd, ymax = Browngreatmed+Browngreatsd), fill = colors[2], alpha = 0.2) + 
   #  # geom_ribbon(aes(ymin = Levygreatmed-Levygreatsd, ymax = Levygreatmed+Levygreatsd), fill = colors[3], alpha = 0.2) + 
   #  # geom_ribbon(aes(ymin = CRWgreatmed-CRWgreatsd, ymax = CRWgreatmed+CRWgreatsd), fill = colors[4], alpha = 0.2) + 
   #  # geom_ribbon(aes(ymin = JDgreatmed-JDgreatsd, ymax = JDgreatmed+JDgreatsd), fill = colors[5], alpha = 0.2) +
   #  geom_point(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szp) + geom_line(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szl) + 
   #  geom_point(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muCRWless), color = colors[5], alpha = 0.6, size = szp) + geom_line(aes(y = muCRWless), color = colors[5], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muJDless), color = colors[6], alpha = 0.6, size = szp) + geom_line(aes(y = muJDless), color = colors[6], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muAllless), color = colors[1], alpha = 0.6, size = szp) + geom_line(aes(y = muAllless), color = colors[1], alpha = 0.6, size = szl) +
   #  theme_bw() + labs(x = "", y = "Fraction of significant tests") + coord_cartesian(clip = 'off')
   #
   #Qbelow <- ggplot(data = ntrackBS, aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), limits = c(1, NA)) +
   #  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
   #  geom_ribbon(aes(ymin = pmax(0,muAllgreat-3*sdAllgreat), ymax = pmin(1,muAllgreat+3*sdAllgreat)), fill = colors[1], alpha = 0.2) +
   #  geom_ribbon(aes(ymin = pmax(0,muBrowngreat-3*sdBrowngreat), ymax = pmin(1,muBrowngreat + 3*sdBrowngreat)), fill = colors[2], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muLevygreat-3*sdLevygreat), ymax = pmin(1,muLevygreat+3*sdLevygreat)), fill = colors[3], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muCRWgreat-3*sdCRWgreat), ymax = pmin(1,muCRWgreat+3*sdCRWgreat)), fill = colors[5], alpha = 0.2) + 
   #  geom_ribbon(aes(ymin = pmax(0,muJDgreat-3*sdJDgreat), ymax = pmin(1,muJDgreat+3*sdJDgreat)), fill = colors[6], alpha = 0.2) +
   #  # geom_ribbon(aes(ymin = Allgreatmed-Allgreatsd, ymax = Allgreatmed+Allgreatsd), fill = colors[1], alpha = 0.2) +
   #  # geom_ribbon(aes(ymin = Browngreatmed-Browngreatsd, ymax = Browngreatmed+Browngreatsd), fill = colors[2], alpha = 0.2) + 
   #  # geom_ribbon(aes(ymin = Levygreatmed-Levygreatsd, ymax = Levygreatmed+Levygreatsd), fill = colors[3], alpha = 0.2) + 
   #  # geom_ribbon(aes(ymin = CRWgreatmed-CRWgreatsd, ymax = CRWgreatmed+CRWgreatsd), fill = colors[4], alpha = 0.2) + 
   # # geom_ribbon(aes(ymin = JDgreatmed-JDgreatsd, ymax = JDgreatmed+JDgreatsd), fill = colors[5], alpha = 0.2) +
   #  geom_point(aes(y = muBrowngreat), color = colors[2], alpha = 0.6, size = szp) + geom_line(aes(y = muBrowngreat), color = colors[2], alpha = 0.6, size = szl) + 
   #  geom_point(aes(y = muLevygreat), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muLevygreat), color = colors[3], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muCRWgreat), color = colors[5], alpha = 0.6, size = szp) + geom_line(aes(y = muCRWgreat), color = colors[5], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muJDgreat), color = colors[6], alpha = 0.6, size = szp) + geom_line(aes(y = muJDgreat), color = colors[6], alpha = 0.6, size = szl) +
   #  geom_point(aes(y = muAllgreat), color = colors[1], alpha = 0.6, size = szp) + geom_line(aes(y = muAllgreat), color = colors[1], alpha = 0.6, size = szl) +
   #  theme_bw() + labs(x = "Number of tracks", y = "Fraction of significant tests") + coord_cartesian(clip = 'off')
    
    
    # pdf(file = "Boosted_pvals_K0_S0.pdf", width = 3, height = 5, pointsize = 16)
    #iff(paste(glue("Boosted_pvals_K{K}_S{S}_less_{CUE}.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
    #rint(Qabove)
    #ev.off()
    #iff(paste(glue("Boosted_pvals_K{K}_S{S}_greater_{CUE}.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
    #rint(Qbelow)
    #ev.off()
    
    Allbootstrapestimates <- rbind(Allbootstrapestimates, ntrackBS)
    Alltesttracks <- rbind(Alltesttracks, testtracks)
    
  }
  
  # save(Allbootstrapestimates, file = "temposave.RData")
#  save(Allbootstrapestimates,file = paste(glue("/lustre/scratch/jpinti/Accuracy/Boostrapestimates_K{K}_{CUE}.RData")))
save(Allbootstrapestimates, Alltesttracks, file = paste(glue("/lustre/scratch/jpinti/Accuracy/Boostrapestimates_{CUE}_K{K}_v4.RData")))
  
}


