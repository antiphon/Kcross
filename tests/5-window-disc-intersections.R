# Test the ball intersection thingy

library(spatstat)

library(devtools)
load_all(".")


x <- runifpoint(10)
r <- 0.1

out <- intersect_bbox_with_disc(x, r)


symbols(x, circles=rep(r, x$n), inches=F, xlim=0:1, ylim=0:1)
rect(bbox[1,1],bbox[1,2],bbox[2,1],bbox[2,2])

rel <- rbind(relArea = out$Area / (pi*r^2), relLEngth=out$Length/(2*pi*r))

print(rel)

# seems to work
