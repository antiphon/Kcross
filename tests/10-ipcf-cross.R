# test local pcf

library(devtools)
load_all(".")
library(spatstat)


x <- rpoispp(100, win = square(2))
r <- seq(0, .3, l=51)
e <- ipcf_local_box(x, r, int = int <- rep(intensity(x), x$n), adjust=1)
f <- localpcf(x, rmax = max(r), nr = 51)
fv <- as.data.frame(f)

kk <- 1/area(x) #epa_integral(r*sqrt(4/x$n)/0.15)


plot(r, colMeans(e)/kk, ylim=c(0,4), "l")
for(i in 1:10)
  lines(r, colMeans(e[(i*40):x$n,])/kk, col=i)


abline(h=1)
#apply(fv[,1:30], 2, lines, x=r, col="gray90")

#apply(e[1:30,], 1, lines, x=r, col="gray60")


