#' Long format from Array
#' @param a Array of the form p x p x m
#' @param r optional scale vector used in computation
#' @details
#' The functions in this package return arrays of dimensions p x p x m, where p = number of types, m = number of scales.
#' This function turns this such arrays into data frames with columns (i, j, r, v). If `r` is provided,
#' the `r` in the output should match. Otherwise, output `r` is the index.
#'
#' @return data.frame
#' @export
array_to_df <- function(a, r) {
  d <- dim(a)
  c(d)
  o <- data.frame(i = rep(1:d[1], length = prod(d[2:3])),
                  j = rep(1:d[2], each = d[1]),  # recycle
                  r = rep(1:d[3], each = prod(d[1:2])),
                  value = c(a))
  if(!missing(r)) o[, 3] <- r[o[, 3]]
  o
}
