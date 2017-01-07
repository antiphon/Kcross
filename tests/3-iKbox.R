# test inhom K

library(devtools)
load_all(".")
library(spatstat)

x <- rpoispp(100)

int <- density(x, at ="points")

k <- Kinhom(x, lambda = int)
k2 <- iK_box(x, r = k$r, int = int)

plot(k$trans, col=2)
lines(k2)


k0 <- K_box(x, r = k$r)
