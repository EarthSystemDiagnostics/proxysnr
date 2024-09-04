#' Check if argument is a spectral object
#'
#' Check passes if object is a named list with the two elements `freq` and
#' `spec` which have equal length, otherwise an error is issued informing about
#' the required object structure.
#' 
#' @param x an object to test.
#' @return \code{TRUE} invisibly if test passes.
#'
#' @author Thomas Münch
#' @noRd
#'
check.if.spectrum <- function(x) {

  nm <- deparse(substitute(x))
  msg <- paste0("`", nm, "` ", "must be a list with elements ",
               "`freq` and `spec` of equal length.")
  
  if (!is.list(x)) stop(msg, call. = FALSE)

  if (!all(utils::hasName(x, c("freq", "spec"))) | stats::sd(lengths(x)) > 0)
    stop(msg, call. = FALSE)

  invisible(TRUE)

}

#' Test for spectral object
#'
#' Function returns \code{TRUE} if object is a named list of the equal-length
#' elements `freq` and `spec`.
#'
#' @param x an object to test whether it is as spectral object.
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @author Thomas Münch
#' @noRd
#'
is.spectrum <- function(x) {

  (tryCatch(check.if.spectrum(x), error = function(cond) {FALSE}))

}

#' Check if `simproxyage` package is installed
#'
#' @param stop.on.false logical whether to stop when package is not installed;
#'   defaults to `FALSE` which only prints the message with the installation
#'   instructions.
#' @return a logical signalling whether the package is installed.
#'
#' @author Thomas Münch
#' @noRd
#'
check.simproxyage <- function(stop.on.false = FALSE) {

  has.simproxyage <- requireNamespace("simproxyage", quietly = TRUE)

  if (has.simproxyage) {

    message("\nChecking simproxyage availability... ok.")

  } else {

    message("\nChecking simproxyage availability...")

    msg <- paste0(
      "Package not found.\n",
      "Install via ",
      "`remotes::install_github(\"EarthSystemDiagnostics/simproxyage\")`\n",
      "or download package repository from ",
      "<https://doi.org/10.5281/zenodo.2025204>\n",
      "and install via `devtools::install()`.\n"
    )

    if (stop.on.false) {
      cat("\n")
      msg <- paste0(msg, "\n")
      stop(msg, call. = FALSE)
    } else {
      msg <- paste0("\n", msg)
      message(msg)
    }
  }

  invisible(has.simproxyage)

}
