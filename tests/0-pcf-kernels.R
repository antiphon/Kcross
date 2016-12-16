# Check box kernel pcf

library(devtools)
load_all(".")

x <- rpoispp(10)


bw <- 0.15/sqrt(10)

r <- seq(0, 0.3, l=50)

ge <- pcf_box(x, r, bw = bw, kern = 0)
gb <- pcf_box(x, r, bw = bw, kern = 1)


plot(ge, ylim=c(0,2))
points(gb, col=2)
abline(h=1)
