# Check K-cross in 3d

library(devtools)
load_all(".")

#
p <- 10
d <- 2
n <- p*d*100
x <- list( x = matrix(runif(d * n, 0, 2), ncol=d), marks = sample(1:p, n, replace=T), bbox = rbind(0, rep(2,d)))
r <- seq(0, 0.25, l = 100)
kest <- K_cross_all_box(x, r = r)

# theo
theo <- pi^(d/2)*r^d / gamma(d/2+1)


plot(r, kest[1,2,], ylim=range(kest) ,type="l")
for(i in 1:p)
  for(j in 1:p) lines(r, kest[i,j,])

lines(r,theo, col=2, lwd=2)
