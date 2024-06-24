#' Test for spectral object
#'
#' Test passes if object is a named list with the two elements `freq` and
#' `spec` which have equal length, otherwise an error is issued.
#' 
#' @param x an object to test.
#' @return \code{TRUE} invisibly if test passes.
#'
#' @author Thomas Münch
#' @noRd
#'
is.spectrum <- function(x) {

  nm <- deparse(substitute(x))
  msg <- paste0("`", nm, "` ", "must be a list with elements ",
               "`freq` and `spec` of equal length.")
  
  if (!is.list(x)) stop(msg, call. = FALSE)

  if (!all(utils::hasName(x, c("freq", "spec"))) | stats::sd(lengths(x)) > 0)
    stop(msg, call. = FALSE)

  invisible(TRUE)

}
