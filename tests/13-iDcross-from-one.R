# test inhom D cross from one type only

library(devtools)
load_all(".")
library(spatstat)
data(lansing)
x <- unique(lansing)
xl <- split(x)

library(Kdirectional)
int1l <- lapply(xl, intensity_at_points, normalise=!TRUE, bw = .2)
int0 <- rep(intensity(x), table(x$marks))

int1 <- unlist(int1l)

int <- int0

# this is the computation
r <- seq(0, 0.2, l = 50)

# compare to spatstat:
m <- levels(x$marks)
mn <- length(m)
mn <- 3
par(mfrow=c(mn,mn))

k3 <- iD_cross_all_box(x, r = r, int = int)


for(i in 1:mn){
  k2 <- iD_cross_all_box_from(x, r = r, int = int, from = i)
  for(j in 1:mn){
  k <- Gcross(x, i = m[i], j = m[j], r = r, correction = "rs")
  plot(r, k$rs, col=2, main=paste(m[i],m[j]), type="l", lwd=2)
  lines(r, k2[,j])
  lines(r, k3[i,j,], col=4, lty=2)
  }
}
