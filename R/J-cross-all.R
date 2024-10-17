# dev J-cross faster

#library(spatstat)
#library(Rcpp)

#sourceCpp("D_cross_all_box.cpp")


D_cross_all_box <- function(x, r) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  intensities <- intensity(x)
  V <- c_D_cross_2d_box(xy, bb, types, intensities, r)
  #
  # compile the pcf's
  G <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    G[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes)
  }
  G
}


# Quick'n'dirty empty space function for a mv pp using spatstat
F_all_box <- function(x, r) {
  xl <- split(x)
  G <- NULL
  G <- sapply(xl, function(x) Fest(x, r=r, correction="km")$km)
  t ( G )
}

# quick cross J, border correction
# assumes rectangular window, uses Fest for denominator
J_cross_all_box <- function(x, r) {
  # compute D cross
  D <- D_cross_all_box(x, r)
  # empty space
  G <- F_all_box(x, r)
  # J
  J <- D
  ntypes <- length(unique(x$marks))
  for(i in 1:ntypes){
    for(j in i:ntypes){
      J[i,j,] <- ( 1 - D[i,j,])/(1 - G[j,])
      J[j,i,] <- ( 1 - D[j,i,])/(1 - G[i,])
    }
  }
  J
}



#############################
# Cross D using spatstat
D_cross_all_sp <- function(x, r) {
  ntypes <- length(unique(x$marks))
  D <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:ntypes){
    for(j in i:ntypes){
      D[i,j,] <- Gcross(x, i, j, r = r, correction = "han")$han
      D[j,i,] <- Gcross(x, j, i, r = r, correction = "han")$han
    }
  }
  D
}

# Cross J using spatstat
J_cross_all_sp <- function(x, r) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  J <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:ntypes){
    for(j in i:ntypes){
      J[i,j,] <- Jcross(x, i, j, r = r, correction = "han")$han
      J[j,i,] <- Jcross(x, j, i, r = r, correction = "han")$han
    }
  }
  J
}

# Check D cross
# if(0){
#   library(spatstat)
#   set.seed(1)
#   N <- 400
#   M <- 4
#   x <- setmarks(rpoint(N), factor(sample(1:M, N, T)))#
#   r <- seq(0, .15, l = 50)
#   #t0 <- system.time( g <- Gcross(x, i=1, r=r, j = 2) )
#   t0 <- system.time( g <- D_cross_all_sp(x, r = r) )
#   t1 <- system.time( gg <- D_cross_all_box(x, r = r) )
#   par(mfrow=c(1,1))
#   plot(r, g[1,2,], "l", ylim =c(0,2), lwd=3)
#   lines(r, g[2,1,], col=1, lty=2, lwd=3)
#   lines(r, gg[1,2,], col=3, lty=2, lwd=2)
#   lines(r, gg[2,1,], col=3, lty=2, lwd=2)
#   abline(h = 1, col=2)
#   print(range(gg))
#   print(rbind(t0,t1))
# }
#
#
# # Check J-cross
# if(0){
#   library(spatstat)
#   set.seed(1)
#   N <- 400
#   M <- 4
#   x <- setmarks(rpoint(N), factor(sample(1:M, N, T)))#
#   r <- seq(0, .1, l = 80)
#   t0 <- system.time( J <- J_cross_all_sp(x, r = r) )
#   t1 <- system.time( JJ <- J_cross_all_box(x, r = r) )
#   par(mfrow=c(1,1))
#   plot(r, J[1,2,], "l", ylim =c(0,2), lwd=3)
#   lines(r, J[2,1,], "l", ylim =c(0,2), lwd=3, col=2, lty=2)
#   lines(r, JJ[1,2,], col=3, lty=2, lwd=2)
#   lines(r, JJ[2,1,], col=4, lty=2, lwd=2)
#   abline(h = 1, col=2)
#   print(rbind(t0,t1))
#   # close enough. Different b-corrections lead to different values.
# }
