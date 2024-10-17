# check quick ISAR

# test inhom D cross
library(devtools)
load_all(".")
library(spatstat)
library(sseg)
library(Kdirectional)

data(lansing)
x <- unique(lansing)
xl <- split(x)

int <- lapply(xl, intensity_at_points, normalise=!TRUE, bw = .2)
int0 <- rep(intensity(x), table(x$marks))


int <- unlist(int)

int1 <- int0

# this is the computation
r <- seq(0, 0.2, l = 50)
k3 <- iISAR_all_box(x, r = r, int = int1)

m <- levels(x$marks)
mn <- length(m)
par(mfrow=c(mn/3,3))

for(i in 1:mn){
  #k <- ISAR_homog(x, i = i,  r = r)#, intensity = int1)
  k <- ISAR_inhom(x, i = i,  r = r, intensity = int1)
  k2 <- iISAR_all_box_from(x, r, i, int1)
  plot(r, k[[2]], col=2, main=paste(m[i]), type="l", ylim=c(0,mn-1), lwd=2)
  lines(r, k2, lty=2)
  lines(r, k3[,i], lty=3, col=4)
}

