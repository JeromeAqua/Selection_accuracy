library("ggplot2")
library("glue")

CUE <- "ftle_GAPS"
K <- 0
# load(glue("../Boostrapestimates_{CUE}_K{K}_v4.RData"))
load(glue("../Boostrapestimates_{CUE}_K{K}.RData"))

##Plot D estimates
Colnames <- c("CUE","kappa", "sigma", "ntracks", "mean_DLess",  "SD_DLess", "n_DLess", 
              "mean_DGreat", "SD_DGreat","n_DGreat") 
Dest <- setNames(data.frame(matrix(ncol = length(Colnames), nrow = 0)), Colnames)
pvalthresh <- 0.05
szl <- 1
szp <- 2


for (s in unique(Alltesttracks$sigma)){
  for (N in unique(Alltesttracks$ntracks)){
mdl <-    mean(Alltesttracks$DallLess[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallLess<pvalthresh]) 
sddl <-     sd(Alltesttracks$DallLess[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallLess<pvalthresh]) 
ndl <-  length(Alltesttracks$DallLess[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallLess<pvalthresh]) 

mdg <-    mean(Alltesttracks$DallGreat[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallGreat<pvalthresh]) 
sddg <-     sd(Alltesttracks$DallGreat[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallGreat<pvalthresh]) 
ndg <-  length(Alltesttracks$DallGreat[Alltesttracks$ntracks==N & Alltesttracks$sigma==s & Alltesttracks$pvalallGreat<pvalthresh]) 

V <- c(CUE, K, s, N, mdl, sddl, ndl, mdg, sddg, ndg)

Dest[nrow(Dest)+1,] <- V
    
  }
}

Dest$ntracks <- as.numeric(Dest$ntracks)
Dest$kappa <- as.numeric(Dest$kappa)
Dest$sigma <- as.numeric(Dest$sigma)

Dest$SD_DLess <- as.numeric(Dest$SD_DLess)
Dest$SD_DGreat <- as.numeric(Dest$SD_DGreat)

Dest$mean_DLess <- as.numeric(Dest$mean_DLess)
Dest$mean_DGreat <- as.numeric(Dest$mean_DGreat)

Dest$n_DLess <- as.numeric(Dest$n_DLess)
Dest$n_DGreat <- as.numeric(Dest$n_DGreat)

colors <- c("black", "gray47", "gray66")

S1 <- 0
S2 <- 50
S3 <- 111
ymax <- .5

Dless <- ggplot(data = Dest[Dest$sigma==S1,], aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), labels = c("","","",""),limits = c(1, NA)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1, .2, .3, .4, .5),labels = c("","","","","",""), limits = c(0,ymax)) +
  geom_ribbon(aes(ymin = pmax(0,mean_DLess-3*SD_DLess), ymax = pmin(ymax,mean_DLess+3*SD_DLess)), fill = colors[1], alpha = 0.2) +
  geom_point(aes(y = mean_DLess), color = colors[1], alpha = 1, size = szp) +geom_line(aes(y = mean_DLess), color = colors[1], alpha = 1, size = szl) +
  geom_ribbon(data = Dest[Dest$sigma==S2,],aes(ymin = pmax(0,mean_DLess-3*SD_DLess), ymax = pmin(ymax,mean_DLess+3*SD_DLess)), fill = colors[2], alpha = 0.2) +
  geom_point(data = Dest[Dest$sigma==S2,],aes(y = mean_DLess), color = colors[2], alpha = 1, size = szp) +geom_line(data = Dest[Dest$sigma==S2,],aes(y = mean_DLess), color = colors[2], alpha = 1, size = szl)+
  geom_ribbon(data = Dest[Dest$sigma==S3,],aes(ymin = pmax(0,mean_DLess-3*SD_DLess), ymax = pmin(ymax,mean_DLess+3*SD_DLess)), fill = colors[3], alpha = 0.2) +
  geom_point(data = Dest[Dest$sigma==S3,],aes(y = mean_DLess), color = colors[3], alpha = 1, size = szp) +geom_line(data = Dest[Dest$sigma==S3,],aes(y = mean_DLess), color = colors[3], alpha = 1, size = szl) + 
  theme_bw() + labs(x = "", y = "") #+ coord_cartesian(clip = 'off')


Dgreat <- ggplot(data = Dest[Dest$sigma==S1,], aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), labels = c("","","",""),limits = c(1, NA)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1, .2, .3, .4, .5), labels = c("","","","","",""), limits = c(0,ymax)) +
  geom_ribbon(aes(ymin = pmax(0,mean_DGreat-3*SD_DGreat), ymax = pmin(ymax,mean_DGreat+3*SD_DGreat)), fill = colors[1], alpha = 0.2) +
  geom_point(aes(y = mean_DGreat), color = colors[1], alpha = 1, size = szp) +geom_line(aes(y = mean_DGreat), color = colors[1], alpha = 1, size = szl) +
  geom_ribbon(data = Dest[Dest$sigma==S2,],aes(ymin = pmax(0,mean_DGreat-3*SD_DGreat), ymax = pmin(ymax,mean_DGreat+3*SD_DGreat)), fill = colors[2], alpha = 0.2) +
  geom_point(data = Dest[Dest$sigma==S2,],aes(y = mean_DGreat), color = colors[2], alpha = 1, size = szp) +geom_line(data = Dest[Dest$sigma==S2,],aes(y = mean_DGreat), color = colors[2], alpha = 1, size = szl)+
  geom_ribbon(data = Dest[Dest$sigma==S3,],aes(ymin = pmax(0,mean_DGreat-3*SD_DGreat), ymax = pmin(ymax,mean_DGreat+3*SD_DGreat)), fill = colors[3], alpha = 0.2) +
  geom_point(data = Dest[Dest$sigma==S3,],aes(y = mean_DGreat), color = colors[3], alpha = 1, size = szp) +geom_line(data = Dest[Dest$sigma==S3,],aes(y = mean_DGreat), color = colors[3], alpha = 1, size = szl) + 
  theme_bw() + labs(x = "", y = "") #+ coord_cartesian(clip = 'off')


tiff(paste(glue("Dless_{CUE}_K{K}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(Dless)
dev.off()
tiff(paste(glue("Dgreat_{CUE}_K{K}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(Dgreat)
dev.off()
