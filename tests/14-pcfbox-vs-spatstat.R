# pcf vs spatstat
library(Kcross)
library(spatstat)

x <- rpoispp(100, nsim=100)
r <- seq(0,.2, l=41)
bw <- 0.015
g <- sapply(lapply(x, pcf, r = r, bw=bw/sqrt(5), divisor="d", kernel="e"), getElement, "trans")
ns <- sapply(x, getElement,"n")
g1 <- sapply(x, pcf_box, r=r,  bw=bw, kernel_correction=1)



par(mfrow=c(2,2))
plot(r, rowMeans(g), ylim=c(0,1.2))
plot(r, rowMeans(g1), ylim=c(0,1.2))

plot(NA, xlim=range(r), ylim=range(-.1,1)*.3)
for(i in 1:ncol(g)) {
  lines(r, g[,i]-g1[,i])
}

