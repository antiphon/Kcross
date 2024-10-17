#' quick cross K, translation border correction
#' assumes rectangnular window
#' @param x point pattern, multivariate
#' @param r range vector
#' @param int vector of intensities per point
#' @param correction use or not use translation correction
#' @param verb print some iterations
#' @export


iK_cross_all_box <- function(x, r, int, correction = TRUE, verb = FALSE) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  m <- "intensities needed. Vector, int per point"
  if(missing(int)) stop(m)
  if(length(int) != x$n) stop(m)
  # go
  V <- c_iK_cross_2d_box(xy, bb, types, ntypes, int, r, correction, verb)
  #
  # compile the K's
  K <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    K[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes)
  }
  # done
  K
}


#' quick K, border correction
#' assumes rectangnular window
#' @param x point pattern
#' @param r range vector
#' @param int vector of intensities per point
#' @param correction use or not use translation correction
#'
#' @export

iK_box <- function(x, r, int, correction = TRUE) {
  m <- "intensities needed. Vector, int per point"
  if(missing(int)) stop(m)
  if(length(int) != x$n) stop(m)
  xy <- as.matrix(coords(x))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  # go
  V <- c_iK_2d_box(xy, bb, int, r, correction)
  #
  S <- area(x$window)/sum(1/int)
  # done
  V*S
}



# Check
if(0){
  library(spatstat)
  W <- square(2)
  set.seed(1)
  n <- 300
  M <- 1
  #x <- setmarks(rpoint(n, win = W), factor(sample(1:M, n, T)))#
  x <- rpoint(n, window = square(2))
  int <- (den <- density(x)/M)[x]
  r <- seq(0, .3, l = 20)
  t0 <- system.time( k <- Kinhom(x, lambda = int, r=r) )
  t1 <- system.time( kk <- iK_box(x, r = r, int = int) )
  par(mfrow=c(1,1))
  o <- pi*r^2
  plot(r, k$trans - o, "l", lwd=3, ylim=c(-1,1)*.01)
  lines(r, kk-o, col=3, lty=2, lwd=2)
  print(rbind(t0,t1))
}

if(0){
  library(spatstat)
  W <- square(2)
  set.seed(14)
  n <- 700
  M <- 2
  x <- setmarks(rpoint(n, win = W), factor(sample(1:M, n, T)))#
  #x <- rpoint(n, window = square(2))
  int <- (den <- density(x)/M)[x]
  ints <- split(int, x$marks)
  i1 <- ints[[1]]
  i2 <- ints[[2]]
  r <- seq(0, .3, l = 20)
  #
  t0 <- system.time( k <- Kcross.inhom(x, i=1, j=2,
                                       lambdaI = i1,
                                       lambdaJ = i2,
                                       r=r) )
  t1 <- system.time( kk <- iK_cross_all_box(x, r = r, int = int) )
  par(mfrow=c(1,1))
  o <- pi*r^2
  plot(r, k$trans - o, "l", lwd=3, ylim=c(-1,1)*.04)
  lines(r, kk[1,2,]-o, col=3, lty=2, lwd=2)
  print(rbind(t0,t1))
}



