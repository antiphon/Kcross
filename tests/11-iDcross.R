# test inhom D cross

library(devtools)
load_all(".")
library(spatstat)
data(lansing)
x <- unique(lansing)
xl <- split(x)

#ints <- lapply(xl, density, at = "points", leaveoneout=F)
# std:
#V <- area(x)
#ints <- lapply(ints, function(v)  v * V/sum(1/v) )
#xl <- x

library(Kdirectional)
int <- lapply(xl, intensity_at_points, normalise=!TRUE, bw = .2)
int0 <- rep(intensity(x), table(x$marks))

int <- unlist(int)
# this is the computation
r <- seq(0, 0.2, l = 50)
k2 <- iD_cross_all_box(x, r = r, int = int)

# compare to spatstat:
m <- levels(x$marks)
mn <- length(m)
par(mfrow=c(mn,mn))

for(i in 1:mn)for(j in 1:mn){
  k <- Gcross(x, i = m[i], j = m[j], r = r, correction = "rs")
  plot(r, k$rs, col=2, main=paste(m[i],m[j]), type="l")
  lines(r, k$theo, col=3)
  lines(r, k2[i,j,])
}
