library("ggplot2")
library("glue")

CUE <- "ftle_GAPS"
K <- 0
# load(glue("../Boostrapestimates_{CUE}_K{K}_v4.RData"))
load(glue("../Boostrapestimates_{CUE}_K{K}.RData"))

S <- 111

colors = c("black","red","blue","chartreuse3")
szl <- 1
szp <- 2

Qabove <- ggplot(data = Allbootstrapestimates[Allbootstrapestimates$sigma==S,], aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), labels = c("","","",""),limits = c(1, NA)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6,0.8,1), limits = c(0,1)) +
  geom_ribbon(aes(ymin = pmax(0,muAllless-3*sdAllless), ymax = pmin(1,muAllless+3*sdAllless)), fill = colors[1], alpha = 0.2) +
  geom_ribbon(aes(ymin = pmax(0,muBrownless-3*sdBrownless), ymax = pmin(1,muBrownless + 3*sdBrownless)), fill = colors[2], alpha = 0.2) + 
  # geom_ribbon(aes(ymin = pmax(0,muLevyless-3*sdLevyless), ymax = pmin(1,muLevyless+3*sdLevyless)), fill = colors[3], alpha = 0.2) + 
  geom_ribbon(aes(ymin = pmax(0,muCRWless-3*sdCRWless), ymax = pmin(1,muCRWless+3*sdCRWless)), fill = colors[3], alpha = 0.2) + 
  geom_ribbon(aes(ymin = pmax(0,muJDless-3*sdJDless), ymax = pmin(1,muJDless+3*sdJDless)), fill = colors[4], alpha = 0.2) +
  geom_point(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szp) + geom_line(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szl) + 
  # geom_point(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szl) +
  geom_point(aes(y = muCRWless), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muCRWless), color = colors[3], alpha = 0.6, size = szl) +
  geom_point(aes(y = muJDless), color = colors[4], alpha = 0.6, size = szp) + geom_line(aes(y = muJDless), color = colors[4], alpha = 0.6, size = szl) +
  geom_point(aes(y = muAllless), color = colors[1], alpha = 0.6, size = szp) + geom_line(aes(y = muAllless), color = colors[1], alpha = 0.6, size = szl) +
  theme_bw() + labs(x = "", y = "Fraction of significant tests") + coord_cartesian(clip = 'off')

Qbelow <- ggplot(data = Allbootstrapestimates[Allbootstrapestimates$sigma==S,], aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), limits = c(1, NA)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8), limits = c(0,1)) +
  geom_ribbon(aes(ymin = pmax(0,muAllgreat-3*sdAllgreat), ymax = pmin(1,muAllgreat+3*sdAllgreat)), fill = colors[1], alpha = 0.2) +
  geom_ribbon(aes(ymin = pmax(0,muBrowngreat-3*sdBrowngreat), ymax = pmin(1,muBrowngreat + 3*sdBrowngreat)), fill = colors[2], alpha = 0.2) + 
  # geom_ribbon(aes(ymin = pmax(0,muLevygreat-3*sdLevygreat), ymax = pmin(1,muLevygreat+3*sdLevygreat)), fill = colors[3], alpha = 0.2) + 
  geom_ribbon(aes(ymin = pmax(0,muCRWgreat-3*sdCRWgreat), ymax = pmin(1,muCRWgreat+3*sdCRWgreat)), fill = colors[3], alpha = 0.2) + 
  geom_ribbon(aes(ymin = pmax(0,muJDgreat-3*sdJDgreat), ymax = pmin(1,muJDgreat+3*sdJDgreat)), fill = colors[4], alpha = 0.2) +
  geom_point(aes(y = muBrowngreat), color = colors[2], alpha = 0.6, size = szp) + geom_line(aes(y = muBrowngreat), color = colors[2], alpha = 0.6, size = szl) + 
  # geom_point(aes(y = muLevygreat), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muLevygreat), color = colors[3], alpha = 0.6, size = szl) +
  geom_point(aes(y = muCRWgreat), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muCRWgreat), color = colors[3], alpha = 0.6, size = szl) +
  geom_point(aes(y = muJDgreat), color = colors[4], alpha = 0.6, size = szp) + geom_line(aes(y = muJDgreat), color = colors[4], alpha = 0.6, size = szl) +
  geom_point(aes(y = muAllgreat), color = colors[1], alpha = 0.6, size = szp) + geom_line(aes(y = muAllgreat), color = colors[1], alpha = 0.6, size = szl) +
  theme_bw() + labs(x = "Number of tracks", y = "Fraction of significant tests") + coord_cartesian(clip = 'off')


# pdf(file = "Boosted_pvals_K0_S0.pdf", width = 3, height = 5, pointsize = 16)
tiff(paste(glue("Boosted_pvals_K{K}_S{S}_less_{CUE}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(Qabove)
dev.off()
tiff(paste(glue("Boosted_pvals_K{K}_S{S}_greater_{CUE}_v4.tiff")), units="cm", width=15, height=10, res=800, compression = 'lzw')
print(Qbelow)
dev.off()

# # Plot D value for significant "less" test
# Dless <- ggplot(data = Alltesttracks[Alltesttracks$sigma==S & Alltesttracks$pvalallLess<0.05,], aes(x = ntracks)) + scale_x_continuous( expand = c(0, 0), breaks = c(20,40,60,80), labels = c("","","",""),limits = c(1, NA)) +
#   scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1, .2, .3, .4, .5), limits = c(0,.5)) +
#   # geom_point(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szp) + geom_line(aes(y = muBrownless), color = colors[2], alpha = 0.6, size = szl) + 
#   # geom_point(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szp) + geom_line(aes(y = muLevyless), color = colors[3], alpha = 0.6, size = szl) +
#   # geom_point(aes(y = muCRWless), color = colors[5], alpha = 0.6, size = szp) + geom_line(aes(y = muCRWless), color = colors[5], alpha = 0.6, size = szl) +
#   # geom_point(aes(y = muJDless), color = colors[6], alpha = 0.6, size = szp) + geom_line(aes(y = muJDless), color = colors[6], alpha = 0.6, size = szl) +
#   geom_point(aes(y = DallLess), color = colors[1], alpha = 0.6, size = szp) + #geom_line(aes(y = DAllless), color = colors[1], alpha = 0.6, size = szl) +
#   theme_bw() + labs(x = "", y = "Test statistic") + coord_cartesian(clip = 'off')

