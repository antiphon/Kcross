# Test the itnensity adpative

library(spatstat)

library(devtools)
load_all(".")


x <- runifpoint(200)
r <- seq(0, 0.2, l=10)

int <- intensity_adapted(x, r)

z<-rho_box(x, r)

g <- z/(pi * int^2)

plot(NA, xlim=range(r), ylim=c(0,2))
apply(g, 2, lines, x=r)
