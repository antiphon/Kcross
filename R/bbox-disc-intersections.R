#' Compute the intersection of a disc and its border with bounding box
#'
#' Page 486 in IPSS08
#'
#' @param x point pattern with rectangular window
#' @param r radius (vector)
#'
#' @export


intersect_bbox_with_disc <- function(x, r) {
  bbox <- cbind(x$window$xrange, x$window$yrange)
  xy <- rbind(as.matrix(coords(x)))
  A <- S <- NULL
  for(ri in r){
    v <- c_bbox_disc_intersection(xy, bbox, ri)
    A <- cbind(A, v[,1])
    S <- cbind(S, v[,2])
  }
  list(r = r, Area = A, Length = S)
}
