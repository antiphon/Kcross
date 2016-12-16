# Dev variance estimator for the Guan KPSS test

library(devtools)
library(spatstat)
load_all(".")

b <- 2

W <- as.owin(c(0,b,0,b))
x <- rMatClust(20, 0.1, 10)
bbox <- cbind(c(0,b), c(0,b))

xy <- cbind(x$x, x$y)

h <- c(0.01, 0.05, 0.1, 0.2)
R <- h

l <- intensity(x)

ds1 <- 2 * c_guan_doublesum_2d_box(xy, bbox, bw = h) - l^2 * pi * h^2 + l
ds2 <- l^2 * ( K_box(x, R) - pi * R^2 ) + l

sigma2 <- cbind(r=R, ds1, ds2)
print(sigma2)
