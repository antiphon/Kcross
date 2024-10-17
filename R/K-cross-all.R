#' Quick cross K, rectangle translation border correction
#'
#' Compute all cross-K's quickly for a multivariate point pattern observed in a 2d or 3d rectangular box.
#'
#' @param x list(x=coords, marks=integer or factor, bbox=bounding box, columnwise ranges)
#' @param r range vector for the function to be estimated at
#' @param correction use translation correction?
#' @param from optional: target type "i" to compute the cross-K_ij
#'
#' @details By default the cross-K is computed for all two-type combinations.
#'
#' @return An array of size (p,p,k) where p is the number of types and k is the length of r.
#' For example, the estimated cross-K(r) of the type pair (i,j) is accessed with  [i,j,].
#'
#' @examples
#' x <- list( x = matrix(runif(3 * 1000, 0, 2), ncol=3), marks = sample(1:5, 1000, replace=T), bbox = cbind(c(0,2),c(0,2),c(0,2)))
#' r <- seq(0, 0.25, l = 100)
#' kest <- K_cross_all_box(x, r = r)
#' plot(r, kest[1,2,])
#'
#' @export
#' @import Rcpp
#' @useDynLib Kcross

K_cross_all_box <- function(x, r, correction = TRUE, from = NULL) {
  x <- check_mpp(x, marked = TRUE)
  #
  xy <- x$x
  types <- as.integer(x$marks)
  type_levels <- x$type_levels
  ntypes <- length(type_levels)
  bb <- x$bbox
  Vol <- prod(apply(bb, 2, diff))
  counts <- as.numeric(table(types))
  #
  intensities <- counts/Vol
  int2 <- outer(intensities, intensities, "*")
  # go
  if(is.null(from))
    V <- c_K_cross_2d_box(xy, bb, types, ntypes, r, correction)
  else {
    fromi <- match(from, type_levels)
    if(length(fromi)==0) stop("can not interpret 'from'")
    V <- c_K_cross_2d_box_just_one_from(xy, bb, types, fromi, ntypes, r, correction)
  }
  #
  if(!correction) int2 = int2 * Vol
  # compile the K's
  G <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    G[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes) / int2
  }
  # done
  G
}


# #' Quick cross K, rectangle translation border correction
# #'
# #'
# #' @export
# #' @import spatstat Rcpp
# #' @useDynLib Kcross
#
# K_cross_all_box_2d <- function(x, r, correction = TRUE, from = NULL) {
#   xy <- as.matrix(coords(x))
#   types <-  as.integer(marks(x))
#   ntypes <- length(unique(types))
#   bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
#   intensities <- intensity(x)
#   counts <- as.numeric(table(types))
#   int2 <- outer(intensities, intensities, "*")
#   # go
#   if(is.null(from))
#     V <- c_K_cross_2d_box(xy, bb, types, ntypes, r, correction)
#   else {
#     fromi <- match(from, levels(marks(x)))
#     if(length(fromi)==0) stop("can not interpret 'from'")
#     V <- c_K_cross_2d_box_just_one_from(xy, bb, types, fromi, ntypes, r, correction)
#   }
#   #
#   if(!correction) int2 = int2 * area(x)
#   # compile the K's
#   G <- array(dim = c(ntypes, ntypes, length(r)))
#   for(i in 1:length(r)){
#     G[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes) / int2
#   }
#   # done
#   G
# }


#' Quick K, rectangle translation border correction
#'
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
K_box <- function(x, r, correction = TRUE) {
  x <- check_mpp(x, marked = FALSE)
  #
  xy <- x$x
  bb <- x$bbox
  Vol <- prod(apply(bb, 2, diff))
  #
  int <- nrow(xy)/Vol
  # go
  V <- c_K_2d_box(xy, bb, r, correction)
  #
  # done
  2*V/int^2
}



#' K cross for all, using spatstat
#'
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
K_cross_sp <- function(x, r, mc.cores = 2) {
  require(parallel)
  sp <- levels(droplevels(x$marks))
  nsp <-  length(sp)
  # upper triangle
  K <- NULL
  for(si in 1:nsp) {
    i <- sp[si]
    k <- parallel::mclapply(si:nsp, function(sj) {
      j <- sp[sj]
      Kcross(x, i, j, r=r, correction="trans")$trans
    }, mc.cores=mc.cores)
    K <- rbind(K, do.call(rbind, k))
  }
  # make an array of matrices
  nr <- length(r)
  stack <- function(K){
    A <- Ai <- array(0, dim=c(nsp, nsp, nr))
    for(i in 1:nr) {
      B <- A[,,i]
      B[upper.tri(B,T)] <- K[,i]
      B[lower.tri(B,T)] <- K[,i]
      diag(B) <- diag(B)/2
      rownames(B) <- colnames(B) <- sp
      A[,,i] <- B
    }
    A
  }
  A <- stack(K)
  A
}



#
# # Check
# if(0){
#   library(spatstat)
#   W <- square(2)
#   #source("load-data.R")
#   #x <- pats[[1]]
#   set.seed(1)
#   n <- 2*3000
#   M <- 10
#   x <- setmarks(rpoint(n, win = W), factor(sample(1:M, n, T)))#
#   r <- seq(0, .3, l = 20)
#   t0 <- system.time( k <- K_cross_sp(x, r=r) )
#   t1 <- system.time( kk <- K_cross_all_box(x, r = r) )
#   par(mfrow=c(1,1))
#   o <- pi*r^2
#   plot(r, k[2,1,]-o, "l", lwd=3, ylim=c(-1,1)*.01)
#   lines(r, kk[2,1,]-o, col=3, lty=2, lwd=2)
#   for(i in 3:10){
#     lines(r, k[i,1,]-o, "l", lwd=3)
#     lines(r, kk[i,1,]-o, col=3, lty=2, lwd=2)
#   }
#
#   print(rbind(t0,t1))
# }
#
#
# # Check univ
# if(0){
#   library(spatstat)
#   set.seed(125)
#   n <- 2*3000
#   x <- rpoint(n)
#   r <- seq(0, .3, l = 20)
#   t0 <- system.time( k <- Kest(x, r=r, correction="trans") )
#   t1 <- system.time( kk <- K_box(x, r = r) )
#   par(mfrow=c(1,1))
#   o <- pi*r^2
#   plot(r, k$trans-o, "l", lwd=3, ylim=c(-1,1)*.01)
#   lines(r, kk-o, col=3, lty=2, lwd=2)
#   print(rbind(t0,t1))
# }
#
#
#
# # Check
# if(0){
#   library(spatstat)
#   set.seed(1)
#   n <- 2*1000
#   M <- 2
#   x <- setmarks(rpoint(n), factor(sample(1:M, n, T)))#
#   r <- seq(0, .3, l = 20)
#   t1 <- system.time( kk <- K_cross_all_box(x, r = r) )
#   par(mfrow=c(1,1))
#   o <- pi*r^2
#   plot(r, kk[2,1,]-o, "l", lwd=3, ylim=c(-1,1)*.01)
#   lines(r, kk[1,1,]-o, "l", lwd=3, ylim=c(-1,1)*.01)
#   print(rbind(t0,t1))
# }
