# test inhom K cross

library(devtools)
load_all(".")
library(spatstat)
data(lansing)
x <- lansing
xl <- split(x)
ints <- lapply(xl, density, at = "points")

xl <- superimpose(a=xl[[1]], b=xl[[2]])

int <- unlist(ints)

k <- Kcross.inhom(x, lambdaI = ints[[1]], lambdaJ = ints[[2]])

k2 <- iK_cross_all_box(x, r = k$r, int = int)

plot(k$trans, col=2)
lines(k$theo, col=3)
lines(k2[1,2,])
