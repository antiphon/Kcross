# dev K weird cross
#
# library(spatstat)
# library(Rcpp)
#
# sourceCpp("K_cross_partial_box.cpp")
#

# assumes rectangnular window
K_cross_partial <- function(x, r) {
  xy <- as.matrix(coords(x))
  types <- as.integer(marks(x))
  ntypes <- length( unique(types) )
  counts <- as.numeric(table(types))
  ints <- intensity(x)
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  int <- intensity(x)
  # go
  V <- c_K_cross_partial_box(xy, types, ntypes, counts, ints, bb, r)
  #
  V
}

# Check
if(0){
  library(spatstat)
  #source("load-data.R")
  #x <- pats[[1]]
  set.seed(1)
  M <- 25
  a <- rMatClust(10, .05, 15, nsim = M)
  x <- superimpose(a)

  r <- seq(0.01, .3, l = 10)
  k <- K_cross_partial(x, r=r)
  i <- im(k)
  plot(i)
}
