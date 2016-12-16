# test guans kpss test

library(devtools)
library(spatstat)
load_all(".")

x1 <- rMatClust(30, 0.1, 20, w = as.owin(c(0,2,0,1)))
x2 <- x1[ runif(x1$n) < exp(-2*x1$x)  ]
x1 <- x1[runif(x1$n) < x2$n/x1$n]

# determine suitable pcf
h <- c(0.05, 0.1, 0.2, .25)

# par(mfrow=c(2,2))
# for(i in 1:4){
#   r <- seq(0,0.3,l=50)
#   g2 <- pcf_box(x2, r, kern=0, bw = h[i])
#   plot(r, g2, ylim=c(0,5))
#   abline(h=1)
# }
#



a <- kpss.test(x1, h = h)
b <- kpss.test(x2, h = h)

print(cbind(x=a$statistic, xrej=a$reject, y=b$statistic, yrej=b$reject))
par(mfrow=c(2,1))
plot(a$S); points(x1)
plot(b$S); points(x2)
