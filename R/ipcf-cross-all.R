#' Quick SOIRS Pcf In a Box
#'
#' @param x multitpye pp
#' @param r range vector
#' @param int vector of intensities, for each point
#' @param adjust bandwidth adjustment, bw(i,j)(x,y) = adjust * 0.15 / sqrt(lambda_j(y))
#' @param correction translation correction on/off
#' @param verb verbosity
#'
#' @export

ipcf_cross_all_box <- function(x, r, int, adjust = 1, correction = TRUE, verb = 0) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  # check int
  m <- "intensities needed. Vector, int per point"
  if(missing(int)) stop(m)
  if(length(int) != x$n) stop(m)
  # go
  V <- c_ipcf_cross_2d_box(xy, bb, ntypes, types, int, r, adjust, correction, verb)
  #
  # Normalisation, to ensure intensities follow Campbell.
  vol <- area(x$window)
  Ss <- sapply(split(int, types), function(x) vol/sum(1/x))
  S <- matrix(Ss, ntypes, ntypes)
  # compile the pcf's
  G <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    G[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes) * S
  }

  # done
  G/2
}

####################################################
#' quick univ inhom pcf, translate border correction
#' assumes rectangnular window
#'
#' @export
ipcf_box <- function(x, r, int, bw, adjust = 1, correction = TRUE){
  m <- "int should be a vector, one value per point."
  if(missing(int)) stop(m)
  if(x$n != length(int)) stop(m)
  xy <- as.matrix(coords(x))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  if(missing(bw)) bw <- 0.15 * adjust / sqrt(mean(int))
  #
  S <- area(x$window)/sum(1/int)
  # go
  V <- c_ipcf_2d_box(xy, bb, int, bw, r, adjust, correction)
  #
  S*V
}


####################################################
#' quick univ local inhom pcf, translate border correction
# assumes rectangnular window
#' @export

ipcf_local_box <- function(x, r, int, bw, adjust = 1, correction = TRUE){
  m <- "int should be a vector, one value per point."
  if(missing(int)) stop(m)
  if(x$n != length(int)) stop(m)
  xy <- as.matrix(coords(x))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  if(missing(bw)) bw <- 0.15 * adjust / sqrt(mean(int))
  #
  #S <- area(x$window)/sum(1/int)
  # go
  V <- c_ipcf_local_2d_box(xy, bb, int, bw, r, adjust, correction)
  #
  V
}


####################################################
# TEST
# univariate ig
#
# if(0){
#   library(spatstat)
#   set.seed(12)
#   x <- rMatClust(100,.05, 10)
#   r <- seq(0, .3, l = 50)
#   int <- density(x)[x]
#   t0 <- system.time( g <- pcfinhom(x, r=r, lambda = int,
#                                    divisor ="d", correction ="trans") )
#   t1 <- system.time( gg <- ipcf_box(x, r = r, int = int) )
#   par(mfrow=c(1,1))
#   plot(r, g$trans, "l", lwd=3)
#   lines(r, gg, col=3, lty=2)
#   abline(h = 1, col=2)
#   print(rbind(t0,t1))
#   # close enough
# }
#
#
#
# # multiv
# if(0){
#   library(spatstat)
#   set.seed(12)
#   n <- 700
#   m <- 10
#   x <- setmarks( rpoint(n), factor( sample(1:m, n, replace=T) )  )
#   r <- seq(0, .4, l = 50)
#
#   int <- (den<-density(x, adjust = 2)/m)[x]
#
#   t0 <- system.time( g <- pcfcross.inhom(x, r=r, i=1, j=2, adjust=2,
#                                          lambdaI = den, lambdaJ=den,
#                                    divisor ="d", correction ="trans") )
#
#   t1 <- system.time( gg <- ipcf_cross_all_box(x, r = r, int = int) )
#   par(mfrow=c(1,1))
#   plot(r, g$trans, "l", lwd=3, ylim= c(0,2))
#   lines(r, gg[1,2,]/2, col=3, lty=2)
#   abline(h = 1, col=2)
#   print(rbind(t0,t1))
#   # close enough
# }
#
