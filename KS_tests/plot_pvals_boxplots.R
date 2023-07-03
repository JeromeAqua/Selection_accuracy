## Plot p values as a function of n
library("ggplot2")
library("glue")

CUE <- "ftle_GAPS"
K <- 0
load(glue("../Boostrapestimates_{CUE}_K{K}.RData"))

S <- 0

colors = c("black","red","blue","darkgoldenrod2","purple","chartreuse3")
szl <- 1
szp <- 2

pless <- ggplot(subset(Alltesttracks, sigma ==S), aes(x = as.factor(ntracks), y = pmax(-20,log10(pvalallLess)))) +            # Applying ggplot function
         geom_boxplot() + theme_bw() + coord_cartesian(ylim = c(-20, 0)) +  
         geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red") + 
         geom_hline(yintercept=log10(0.01), linetype="dashed", color = "red")

pgreat <- ggplot(subset(Alltesttracks, sigma ==S), aes(x = as.factor(ntracks), y = pmax(-20,log10(pvalallGreat)))) +            # Applying ggplot function
          geom_boxplot() + theme_bw() +  coord_cartesian(ylim = c(-20, 0)) + 
          geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red") + 
          geom_hline(yintercept=log10(0.01), linetype="dashed", color = "red")


tiff(paste(glue("pvalLess_{CUE}_K{K}_S{S}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(pless)
dev.off()

tiff(paste(glue("pvalGreat_{CUE}_K{K}_S{S}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(pgreat)
dev.off()


