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

  if (!all(utils::hasName(x, c("freq", "spec")))) stop(msg, call. = FALSE)

  if (stats::sd(lengths(x[c("freq", "spec")])) > 0) stop(msg, call. = FALSE)

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

#' Check frequency axis overlap
#'
#' Check if the frequency axis of a given spectral object overlaps with the
#' frequency axis of a target spectrum. This check is useful before
#' interpolating a spectrum onto some target frequency axis.
#'
#' @param x a spectral object.
#' @param target as \code{x} for the target spectrum.
#' @return \code{TRUE} when the target frequency axis falls within the range of
#'   the frequency axis of \code{x}, \code{FALSE} otherwise.
#'
#' @author Thomas Münch
#' @noRd
#'
has.common.freq <- function(x, target) {

  ftarget <- target$freq
  fx <- x$freq

  min(fx) %<=% min(ftarget) & max(fx) %>=% max(ftarget)

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

#' Spectral object window
#'
#' Extract a subset of a given spectral object observed between the specified
#' lower and upper frequencies.
#'
#' @param x a spectral object.
#' @param f.start the lower frequency of interest; the default \code{NULL} uses
#'   the lowest frequency of \code{x}.
#' @param f.end the upper frequency of interest; the default \code{NULL} uses
#'   the uppermost frequency of \code{x}.
#' @return a subset of the spectral object as observed between \code{f.start}
#'   and \code{f.end}. If both \code{f.start} and \code{f.end} are \code{NULL}
#'   the input \code{x} is returned.
#'
#' @author Thomas Münch
#' @noRd
#'
fwindow <- function(x, f.start = NULL, f.end = NULL) {

  if (is.null(f.start) & is.null(f.end)) return(x)

  n <- length(x$freq)
  i.lo <- if (is.null(f.start)) 1 else min(which(x$freq %>=% f.start))
  i.hi <- if (is.null(f.end)) n else max(which(x$freq %<=% f.end))

  x$freq <- x$freq[i.lo : i.hi]
  x$spec <- x$spec[i.lo : i.hi]

  return(x)

}

#' Linear regression on spectral object
#'
#' Perform a linear regression of a spectral object in logarithmic space to fit
#' a power law of the form `alpha * f^{-beta}` where `f` denotes frequency.
#'
#' @param s a spectral object.
#' @param f.start the lower end of the frequency range on which the fit is
#'   estimated; the default \code{NULL} uses the lowest frequency of \code{s}.
#' @param f.end the upper end of the frequency range on which the fit is
#'   estimated; the default \code{NULL} uses the uppermost frequency of
#'   \code{s}.
#' @return a list with the `alpha` and `beta` parameters of the power law.
#'
#' @author Thomas Münch
#' @noRd
#'
fit.powerlaw <- function(s, f.start = NULL, f.end = NULL) {

  s <- fwindow(s, f.start, f.end)

  x <- s$freq
  y <- s$spec

  y[y < 0] <- NA

  coefs <- stats::lm(log(y) ~ log(x)) %>%
    stats::coef()

  list(alpha = exp(coefs[1]), beta = -1 * coefs[2])

}
