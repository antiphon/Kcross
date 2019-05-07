#' Space-Time SOIRS Kcross-PCF
#'
#'
#' @param x point pattern as an  n x (d+2) -matrix of coordinates
#' @param int intensity at points
#' @param r vector of spatial lags
#' @param t vector of temporal lags
#' @param bbox bounding box for translation correction
#' @param sigmas vector c(ss, st) with space and time bandwidth (sd of gaussian), respectively
#' @param do_correction translation correction? only sensible for cuboidal regions
#' @details
#' Input matrix x of dimensions n x (d+2): columns 1:d taken as space dimensions, d+1 taken as the time dimension and d+2 taken as the type.
#'
#' Note that the bandwidths are fixed for all pairs of types.
#' This is not optimal for highly imbalanced patterns.
#'
#' @return array of dimensions (ntype, ntype, nt, nr)
#'
#' @export
ipcf_st_cross_all_box <- function(x, int, r, t, bbox, sigmas, do_correction = TRUE) {
  xy <- as.matrix(x)
  n <- nrow(xy)
  d <- ncol(xy)
  types <-  as.integer(xy[,d])
  ntypes <- length(unique(types))
  # check int
  m <- "intensities needed. Vector, int per point."
  if(missing(int)) stop(m)
  if(length(int) != n) stop(m)
  # check bbox for trans
  if(missing(bbox)) bbox <- apply(xy[,1:(d-1)], 2, range) # in case
  # check bandwidths
  m <- "smoothing bandwidths (gaussian sd's) needed. Vector of length 2."
  if(missing(sigmas)) stop(m)
  if(length(sigmas) != 2) stop(m)
  # go
  V <- c_ipcf_st_cross_2d_box(xy, bbox, ntypes, types,
                              int,
                              r, t,
                              sigmas,
                              as.integer( do_correction) )
  #
  vol <- prod( apply(bbox, 2, diff) )
  nt <- length(t)
  nr <- length(r)
  # scaling
  Ss <- sapply(split(int, types), function(x) vol/sum(1/x))
  S <- matrix(Ss, ntypes, ntypes)
  if(!do_correction) S <- S / vol
  # compile the pcf's
  G <- array(dim = c(ntypes, ntypes, nt, nr))
  for(ir in 1:nr)
    for(it in 1:nt){
      G[,,it,ir] <- matrix(V[[(ir-1)*nt + it]], ncol = ntypes, nrow = ntypes) * S
  }
  #browser()

  # done
  G
}
