#' Calculate signal and noise spectra
#' 
#' Calculate the raw signal and noise spectra, and the corresponding
#' signal-to-noise ratio, as estimated from a core array of \code{n} proxy
#' records, and correct these, where applicable, for the effects of time
#' uncertainty and water vapour diffusion (relevant for firn and ice cores).
#'
#' This function is an implementation of Eq. (4) in Münch and Laepple
#' (2018). While the diffusion correction via a diffusion transfer function is
#' relevant only for diffusing proxies such as stable isotopes from firn and ice
#' cores, this function can be applied to a large set of proxy data since only
#' one of the two transfer functions, or none, need to be supplied; thus, e.g.,
#' it can also be used for non-diffusing proxy data where only time uncertainty
#' is relevant, or for estimating raw signal and noise spectra by supplying no
#' transfer functions at all.
#'
#' @param spectra a list of the raw spectral estimates from a proxy core
#'   array as output from \code{\link{ObtainArraySpectra}}, but sufficient is a
#'   named list of two components giving the \code{mean} and the \code{stack}
#'   spectrum.
#' @param neff the effective number of records (e.g. to account for an expected
#'   spatial correlation of the local noise). Per default set to element
#'   \code{N} in \code{spectra}, otherwise supply it explicitly here.
#' @param diffusion a spectral object (= a list of the equal-length vectors
#'   `freq` and `spec`) of a diffusion transfer function (see also
#'   \code{\link{CalculateDiffusionTF}}) to correct for the effect of diffusion;
#'   i.e., internally the inverse of the transfer function values are used to
#'   correct the spectra (see Eq. 4 in Münch and Laepple, 2018). The default
#'   `NULL` is to apply no correction.
#' @param time.uncertainty as \code{diffusion} a transfer function of the
#'   effect of time uncertainty (see also
#'   \code{\link{CalculateTimeUncertaintyTF}}).
#'
#' @return A list of three components, each of class \code{"spec"}:
#'   \describe{
#'   \item{\code{signal}:}{the raw or corrected signal spectrum;}
#'   \item{\code{noise}:}{the raw or corrected noise spectrum;}
#'   \item{\code{snr}:}{the signal-to-noise ratio as calculated from the previous
#'     components.}
#' }
#'
#' @author Thomas Münch
#' @seealso \code{\link{ObtainArraySpectra}}, \code{\link{CalculateDiffusionTF}},
#'   \code{\link{CalculateTimeUncertaintyTF}}
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @export
#'
SeparateSignalFromNoise <- function(spectra, neff = spectra$N,
                                    diffusion = NULL, time.uncertainty = NULL) {

  # error checking

  if (!is.list(spectra)) stop("`spectra` must be a list.", call. = FALSE)

  if (!all(utils::hasName(spectra, c("mean", "stack"))))
    stop("`spectra` must have elements `mean` and `stack`.", call. = FALSE)

  check.if.spectrum(spectra$mean)
  check.if.spectrum(spectra$stack)

  if (length(spectra$mean$freq) != length(spectra$stack$freq)) {
    stop("`mean` and `stack` must have the same number of spectral estimates.",
         call. = FALSE)
  }

  if (!all(spectra$mean$freq == spectra$stack$freq)) {
    stop("Frequency axes of `mean` and `stack` do not match.")
  }

  if (is.null(neff)) {
    stop("Supply (effective) number of records.")
  }

  if (is.null(diffusion)) {

    dtf.corr <- 1

  } else {

    check.if.spectrum(diffusion)

    if (length(diffusion$freq) != length(spectra$mean$freq)) {
      stop("Length of diffusion correction ",
           "does not match length of spectral estimates.")

    }

    dtf.corr <- 1 / diffusion$spec

  }

  if (is.null(time.uncertainty)) {

    ttf.corr <- 1

  } else {

    check.if.spectrum(time.uncertainty)

    if (length(time.uncertainty$freq) != length(spectra$mean$freq)) {
      stop("Length of time uncertainty correction ",
           "does not match length of spectral estimates.")

    }

    ttf.corr <- 1 / time.uncertainty$spec

  }

  # calculate the signal, noise, and snr (signal/noise) spectra

  N     <- neff
  eta   <- N / (N - ttf.corr)
  mean  <- spectra$mean$spec
  stack <- spectra$stack$spec

  signal <- noise  <- snr <- list()

  signal$freq <- noise$freq <- snr$freq <- spectra$mean$freq

  signal$spec <- eta * ttf.corr * dtf.corr * (stack - mean / N)
  noise$spec  <- eta * dtf.corr * (mean - ttf.corr * stack)

  snr$spec <- signal$spec / noise$spec
  
  # organize output

  class(signal) <- class(noise) <- class(snr) <- "spec"
  
  res <- list(
    signal  = signal,
    noise   = noise,
    snr     = snr)

  return(res)

}
