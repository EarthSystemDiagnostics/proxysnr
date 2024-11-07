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

#' Interpolate spectrum
#'
#' Interpolate a spectrum onto a given frequency axis.
#'
#' @param x a spectral object which is to be interpolated onto the frequency
#'   axis given by the \code{target} spectrum.
#' @param target as \code{x}, used to supply the target frequency axis for the
#'   interpolation.
#' @param num.prec number of decimal places to round the frequency axes in order
#'   to prevent erroneous NA values in the interpolation from floating point
#'   machine representation accuracy. Be careful when changing this value!
#' @return a spectral object with the spectrum in \code{x} interpolated onto the
#'   target frequency axis.
#'
#' @author Thomas Münch
#' @noRd
#'
InterpolateSpectrum <- function(x, target, num.prec = 8) {

  if (!has.common.freq(x, target))
    warning("NAs produced in interpolation as frequency axes do not overlap.",
            call. = FALSE)

  result <- list(
    freq = target$freq,
    spec = stats::approx(round(x$freq, num.prec), x$spec,
                         round(target$freq, num.prec))$y
  )

  class(result) <- "spec"

  return(result)

}

has.array.attribute <- function(x) {

  n <- sprintf(" `%s`.", deparse(substitute(x)))
  a <- attributes(x)

  if (length(a) == 0)
    stop("Attribute `array.par` missing from input object", n, call. = FALSE)

  if (length(a$array.par) == 0)
    stop("Attribute `array.par` missing from input object", n, call. = FALSE)

  if (!all(utils::hasName(a$array.par, c("nc", "nt", "res"))))
    stop("Attribute `array.par` must a named vector with elements ",
         "`nc`, `nt`, `res`.", call. = FALSE)

  if (!checkmate::testNumber(a$array.par[["nc"]], lower = 2, finite = TRUE))
    stop("Element `nc` of `array.par attribute ",
         "(number of proxy records) ",
         "must be a single integer >= 2.", call. = FALSE)

  if (!checkmate::testNumber(a$array.par[["nt"]], lower = 9, finite = TRUE))
    stop("Element `nt` of `array.par attribute ",
         "(number of observations per proxy record) ",
         "must be a single integer > 8.", call. = FALSE)

  if (!checkmate::testNumber(a$array.par[["res"]], lower = 1e-15, finite = TRUE))
    stop("Element `res` of `array.par attribute ",
         "(resolution of proxy records) ",
         "must be a single integer > 0.", call. = FALSE)

}

