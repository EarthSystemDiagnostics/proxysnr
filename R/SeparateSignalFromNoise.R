#' Calculate signal and noise spectra
#' 
#' Calculate the raw signal and noise spectra and the corresponding
#' signal-to-noise ratio spectrum from the spectral estimates of a core array of
#' \code{n} proxy records. Where applicable, these raw results can be corrected
#' for the effects of time uncertainty and diffusion-like smoothing.
#'
#' This function is an implementation of Eq. (4) in Münch and Laepple
#' (2018). While the diffusion transfer function there specifically refers to
#' the diffusional smoothing of stable isotopes from firn and ice cores, it can
#' be interpreted in a much more general sense as a transfer function that
#' describes any smoothing process affecting a proxy record (e.g., bioturbation
#' in marine sediment or biological memory in tree ring records). Therefore,
#' this function can be applied to a large set of proxy data, also because the
#' application of the transfer functions is flexible: e.g., it can be applied on
#' proxy data where only time uncertainty is relevant, or for estimating raw
#' signal and noise spectra by supplying no transfer functions at all.
#'
#' @param spectra a list of the spectral estimates from a proxy core array as
#'   output from \code{\link{ObtainArraySpectra}}, or, as minimum requirement, a
#'   named list (components `mean` and `stack`) supplying the mean spectrum of
#'   the `n` proxy records and the spectrum of the record stacked across the `n`
#'   records.
#' @param neff the effective number of records (e.g. to account for an expected
#'   spatial correlation of the local noise). Per default set to element
#'   \code{N} in \code{spectra}, otherwise supply it explicitly here.
#' @param diffusion a spectral object (= a list of the equal-length vectors
#'   `freq` and `spec`) of a transfer function desribing a diffusion-like proxy
#'   smoothing process (see Details), e.g. diffusion in ice cores (see also
#'   \code{\link{CalculateDiffusionTF}}). Internally, the inverse of the
#'   transfer function values are applied to correct for the smoothing effect on
#'   the estimated signal and noise spectra (see Eq. 4 in Münch and Laepple,
#'   2018). The default `NULL` is to apply no correction.
#' @param time.uncertainty as \code{diffusion} but for a transfer function
#'   that describes the effect of time uncertainty (see also
#'   \code{\link{CalculateTimeUncertaintyTF}} for calculating transfer functions
#'   in the case of layer-counted proxy chronologies) and which is used to
#'   correct the effect it has on the estimated signal spectrum. The default
#'   `NULL` is to apply no correction.
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
