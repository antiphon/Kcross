#' check pattern format
#'
#' @param x pattern candidate
#' @param marked whether to parse the marks as well?
#'
#' @import spatstat
#' @export

check_mpp <- function(x, marked = TRUE){
  if("ppp"%in%is(x)){
    window <- x$window
    if(marked) {
      m <- marks(x)
      if(is.null(m)) stop("marks not found.")
    } else m <- NULL
    x <- list(x=as.matrix(coords(x)), marks = m)
    x$bbox <- cbind(window$xrange, window$yrange)
  }
  merror <- "x should be a list(x = (coordinates), marks = (factor vector), bbox = cbind(x-range, y-range))"
  if(!is(x,"list")){
    if(marked) stop(merror)
    if(is(x, "matrix")) x <- list(x=x, bbox = apply(x, 2, range))
    else stop(paste("Can not parse x of type", is(x) )  )
  }
  else if(is.null(x$bbox)) if(marked) stop(merror) else stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  x$x <- as.matrix(x$x)
  x$bbox <- as.matrix(x$bbox)
  if(ncol(x$x)!=ncol(x$bbox)) stop("x$x and x$bbox have different dimensions.")
  if(marked) {
    if(length(x$marks)  != nrow(x$x)) stop("length(marks) != n(points)")
    x$marks <- factor(x$marks)
    x$type_levels <- levels(x$marks)
  }
  x
}

