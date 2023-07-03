
library(circular)
library(ggplot2)


##Make pretty figure of von mises distribution
x <- seq(-pi,pi,.01)
a <- as.numeric(dvonmises(x,0,2))
b <- as.numeric(dvonmises(x,0,.75))
c <- as.numeric(dvonmises(x,0,.25))
d <- as.numeric(dvonmises(x,0,0))
# e <- as.numeric(dvonmises(x,0,20))

DF <- data.frame(x,a,b,c,d)

ggplot(DF, aes(x = x, y = a)) + geom_line() + geom_line(aes(y=b)) +  geom_line(aes(y=c)) +  geom_line(aes(y=d)) + # geom_line(aes(y=e)) +
  scale_x_continuous( expand = c(0, 0), breaks = c(-pi, -pi/2, 0, pi/2, pi), labels = c("-180","-90","0","90", "180")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,.55)) +
  theme_bw() + labs(x = "", y = "") + coord_cartesian(clip = 'off')



##Make pretty figure of the step length distribution
x <- seq(-80, 80, 0.1)
a <- dnorm(x,0,10/(2+0.5))
b <- dnorm(x,0,10/(.75+0.5))
c <- dnorm(x,0,10/(.25+0.5))
d <- dnorm(x,0,10/(0+0.5))
DF <- data.frame(x,a,b,c,d)

ggplot(DF, aes(x = x, y = a)) + geom_line() + geom_line(aes(y=b)) +  geom_line(aes(y=c)) +  geom_line(aes(y=d)) + # geom_line(aes(y=e)) +
  scale_x_continuous( expand = c(0, 0), breaks = c(-80,-40,0,40,80), labels = c("-80","-40","0","40", "80")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,.101)) +
  theme_bw() + labs(x = "", y = "") + coord_cartesian(clip = 'off')

