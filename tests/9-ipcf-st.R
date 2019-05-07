# Test ipcf_st

devtools::load_all()

# Should work with just spatial as well:
if(0) {
  set.seed(1)
  # test data
  p <- 3
  n <- 300 * p
  L <- 30
  d <- 2
  V <- L^2
  x <- cbind(x=matrix(runif(n*d, 0, L), nc=d),
             t=1, m=sample(1:p, n, replace=TRUE))
  bb <- cbind(0:1 * L, 0:1 * L, 0:1 * 1)
  int <- table(x[,4])[x[,4]]/V
  intt <- n
  intx <- n/L^2
  t <- 0
  r <- seq(0, L*0.5)
  sigmas <- c(.15/sqrt(intx)/5, 1)# / sqrt(5)

  o  <- ipcf_st_cross_all_box(x, int, r, t, bbox = bb, sigmas = sigmas)
  o0 <- pcf_cross_all_box(ppp(x[,1],x[,2],marks =factor(x[,4]),window=square(L)), r=r)

  plot(r, o[1,1,1,], ylim = c(0,5))
  lines(r , o0[1,1,])
}

####
# True spatio-temporal
if(1) {
  set.seed(1)
# test data
  p <- 4
  n <- 300 * p
  L <- 30
  T <- 10
  d <- 2
  V <- L^2 * T
  x <- cbind(x=matrix(runif(n*d, 0, L), nc=d),
             t=runif(n)*T, m=sample(1:p, n, replace=TRUE))
  int <- table(x[,d+2])[x[,d+2]]/V
  intt <- n/T
  intx <- n/L^2/T
  t <- 0
  r <- seq(0, L*0.5)[-1]
  sigmas <- c(.15/sqrt(intx), 2/intt^(1/3))
  #sigmas <- c(.15/sqrt(intx), 1)
  #x[,3] <- 1
  #bb <- cbind(0:1 * L, 0:1 * L, 0:1 * 1)

  o <- ipcf_st_cross_all_box(x, int, r, t, bbox = bb, sigmas = sigmas)

  plot(r, o[1,1,1,], ylim = c(0,5))
  print(o[1,1,1,])

}
