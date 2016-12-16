# test K from_one

library(devtools)
load_all(".")
library(spatstat)

x0 <- x <- rmpoispp(rep(100, 20)*2)
r <- seq(0, .3, l=20)

#Rprof()
t0 <- system.time(  k <- K_cross_all_box(x, r)  )
t1 <- system.time(  k1 <- K_cross_all_box_2d(x, r, from=1)  )
#summaryRprof()

print( all.equal(k[1,,], k1[1,,]) )

print(rbind(t0,t1))
