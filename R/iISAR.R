#' Quick ISAR function, all
#'
#' ISAR, inhomog.
#'
#' @param x point pattern, multivariate
#' @param r range vector
#' @param int vector of intensities per point
#' @param correction use or not use reduced sample border correction
#'
#' @export

iISAR_all_box <- function(x, r, int, correction=TRUE) {
  D <- iD_cross_all_box(x, r, int, correction)
  t(apply(D,3, rowSums) - apply(D, 3, diag))
}

#' Quick ISAR function, target type
#'
#' ISAR, inhomog.
#'
#' @param x point pattern, multivariate
#' @param r range vector
#' @param from target type, integer
#' @param int vector of intensities per point
#' @param correction use or not use reduced sample border correction
#'
#' @export
iISAR_all_box_from <- function(x, r, from, int, correction=TRUE) {
  D <- iD_cross_all_box_from(x, r, int, correction, from = from)
  #browser()
  rowSums(D[,-from])
}


