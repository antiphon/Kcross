# test guans kpss test K version
#
# No difference, so keep the K-version.
#
library(devtools)
library(spatstat)
load_all(".")

x1 <- rMatClust(50, 0.1, 50, w = as.owin(c(0,20,0,10)))
x2 <- x1[ runif(x1$n) < exp(-2*x1$x)  ]
x1 <- x1[runif(x1$n) < x2$n/x1$n]

# determine suitable pcf
h <- c(0.05, 0.1, 0.2, .25)

t0 <- system.time( a <- kpss.test(x2, h = h, useK=F) )
t1 <- system.time( b <- kpss.test(x2, h = h, useK=T) )

print(cbind(x=a$statistic, xrej=a$reject, y=b$statistic, yrej=b$reject))

print(rbind(t0,t1))
