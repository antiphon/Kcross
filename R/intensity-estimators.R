#' Intensity estimator for 2D pattern in rectangle
#'
#' Adapted homogeneous intensity estimators for use with other estimators.
#'
#' @param x point pattern with rectangluar window
#' @param r range parameter in case needed
#' @param type estimator to use.
#'
#' @details type:
#'
#' * "n": none
#' * "b": minus sampling
#' * "v": relative area
#' * "s": relative border length
#'
#' as in IPSS'08 p 193-194.
#' @export

intensity_adapted <- function(x, r, type = "all") {
  n <- x$n
  W <- x$window
  a <- diff(W$xrange)
  b <- diff(W$yrange)
  if(type == "all") type <- c("n","b","v","s")

  gammahat <- function(t) {
    a * b - 2 * t * (a+b)/pi + t^2 /pi
  }
  out <- NULL
  # naive
  if("n" %in% type){
    out <- cbind(out, naive=rep(n/area(W), length(r)))
  }
  # minus sampling
  if("b" %in% type) {
    out <- cbind(out, border=sapply(r, function(ri)  {wr <- erosion(W, ri); x[wr]$n/area(wr)}  ))
  }
  if("v" %in% type | "s" %in% type) {
    V <- intersect_bbox_with_disc(x, r)
  }
  if("v"%in% type){
    A <- V$Area
    v <- colSums(A)
    l <- 0.5 * a*b*r^2 - 2*(a+b)*r^3/(3*pi) + r^4/(4*pi)
    w <- v / (2 * pi * l)
    out <- cbind(out, "rel_area" = w)
  }
  if("s" %in% type) {
    L <- V$Length
    v <- colSums(L) / (2 * pi * r * gammahat(r))
    out <- cbind(out, "rel_length" = v)
  }
  out
}
