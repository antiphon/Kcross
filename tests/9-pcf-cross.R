# test inhom K cross

library(devtools)
load_all(".")
library(spatstat)
data(lansing)
x <- lansing


g <- pcf_cross_all_box(x, seq(0, .3, l = 50), verb=1)

plot(g[1,1,])
lines(g[1,2,])
lines(g[2,2,],col=3)

