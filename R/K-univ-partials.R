# dev K weird

# library(spatstat)
# library(Rcpp)
#
# sourceCpp("K_univ_partial_box.cpp")


# assumes rectangnular window
K_partial <- function(x, r) {
  xy <- as.matrix(coords(x))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  int <- intensity(x)
  # go
  V <- c_K_partial_box(xy, bb, r)
  #
  V/int/x$n
}

# Check
if(0){
  library(spatstat)
  #source("load-data.R")
  #x <- pats[[1]]
  set.seed(1)
  n <- 2*20
  x <-rMatClust(10, .05, 15)
  r <- seq(0, .3, l = 50)
  k <- K_partial(x, r=r)

  par(mfrow=c(2,1))
  plot(r, diag(k), "l")
  lines(r, Kest(x, r=r)$trans-pi*r^2, col=2)
  plot(im(k))
}
