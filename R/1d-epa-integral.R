#' integral of 1d epa kernel
#'
#' bw = 1
#'
#' @param u upper limit of integration
#'
#' @export
#'
epa_integral <- function(u) {
  v <- 0.75 * (u - u*u*u/3.0 + 2.0/3.0)
  v[u > 1] <- 1
  v[u < -1] <- 0
  v
}
