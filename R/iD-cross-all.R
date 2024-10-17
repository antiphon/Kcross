#' quick inhomogeneous cross D, minus border correction
#'
#' Nearest neighbour distance distribution function, assumes rectangnular window.
#'
#' @param x point pattern, multivariate
#' @param r range vector
#' @param int vector of intensities per point
#' @param correction use or not use reduced sample border correction
#'
#' @export


iD_cross_all_box <- function(x, r, int, correction = TRUE) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  m <- "intensities needed. Vector, int per point"
  if(missing(int)) stop(m)
  if(length(int) != x$n) stop(m)
  # Need typewise min intensities
  min_int <- sapply(split(int, types), min)
  # need border distances
  bdist <- bdist.points(x)
  # go
  V <- c_iD_cross_2d_box(xy, bb, types, ntypes, int, min_int, bdist, r, correction)
  #
  # compile the D's
  #browser()
  D <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    D[,,i] <- 1.0 - matrix(V[[i]], ncol = ntypes, nrow = ntypes)
  }
  # done
  D
}

#' quick inhomogeneous cross D, minus border correction, one-to-many
#'
#' Nearest neighbour distance distribution function, assumes rectangnular window.
#'
#' @param x point pattern, multivariate
#' @param r range vector
#' @param from target type, integer in mark-levels
#' @param int vector of intensities per point
#' @param correction use or not use reduced sample border correction
#'
#' @export


iD_cross_all_box_from <- function(x, r, from, int, correction = TRUE) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  m <- "intensities needed. Vector, int per point"
  if(missing(int)) stop(m)
  if(length(int) != x$n) stop(m)
  if(!from %in% types) stop("from not in types")
  # Need typewise min intensities
  min_int <- sapply(split(int, types), min)
  # need border distances
  bdist <- bdist.points(x)
  # go
  V <- c_iD_cross_2d_box_from(xy, bb, types, ntypes, int, min_int, bdist, r, correction, from)
  #
  # compile the D's
  D <- 1 - do.call(rbind, V)
  colnames(D) <- levels(marks(x))
  attr(D, "from") <- levels(marks(x))[from]
  # done
  D
}

