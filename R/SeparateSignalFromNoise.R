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
#' @param spectra a list of the spectral estimates from a proxy core array in
#'   the format as output from \code{\link{ObtainArraySpectra}}.
#' @param neff the effective number of records (`neff` <= `n`, e.g. to account
#'   for an expected spatial correlation of the local noise). Per default
#'   extracted from the "array.par" attribute of `spectra` (see the "Value"
#'   section from \code{\link{?ObtainArraySpectra}}), but you can supply the
#'   `neff` explicitly here to overwrite the default value.
#' @param measurement.noise a measurement noise level for correcting the proxy
#'   noise spectrum: either a single value or a spectral object. In the former
#'   case, the measurement noise is assumed to exhibit a white spectrum and the
#'   given value is interpreted as its total variance; the latter case can be
#'   applied if the measurement noise has a known spectral shape different from
#'   white noise (but the frequency range of the given measurement noise
#'   spectral object must then cover the frequency range of the proxy
#'   spectra). The default `NULL` assumes no measurement noise.
#' @param diffusion a spectral object of a transfer function desribing a
#'   diffusion-like proxy smoothing process (see Details), e.g. diffusion in ice
#'   cores (see also \code{\link{CalculateDiffusionTF}}). Internally, the
#'   inverse of the transfer function values are applied to correct for the
#'   smoothing effect on the estimated signal and noise spectra (see Eq. 4 in
#'   Münch and Laepple, 2018). The default `NULL` is to apply no correction.
#' @param time.uncertainty as \code{diffusion} but for a transfer function
#'   that describes the effect of time uncertainty (see also
#'   \code{\link{CalculateTimeUncertaintyTF}} for calculating transfer functions
#'   in the case of layer-counted proxy chronologies) and which is used to
#'   correct the effect it has on the estimated signal spectrum. The default
#'   `NULL` is to apply no correction.
#'
#' @return A list of three spectral objects:
#'   \describe{
#'   \item{\code{signal}:}{the raw or corrected signal spectrum;}
#'   \item{\code{noise}:}{the raw or corrected noise spectrum;}
#'   \item{\code{snr}:}{the signal-to-noise ratio as calculated from the previous
#'     components;}
#' }
#' with the attribute "array.par": a named vector with information on the proxy
#' record array: number of (effective) records ("nc" = \code{neff}), number of
#' observation points per record ("nt"), and sampling resolution ("res").
#'
#' @author Thomas Münch
#' @seealso \code{\link{ObtainArraySpectra}}, \code{\link{CalculateDiffusionTF}},
#'   \code{\link{CalculateTimeUncertaintyTF}}, `?spec.object` for the definition
#'   of a "proxysnr" spectral object.
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @export
#'
SeparateSignalFromNoise <- function(spectra, neff = NULL,
                                    measurement.noise = NULL,
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

  has.array.attribute(spectra)

  if (is.null(neff)) {

    neff <- attr(spectra, "array.par")[["nc"]]

  } else {

    if (!checkmate::testNumber(neff, lower = 2, finite = TRUE))
      stop("Manually supplied `neff` must be a single integer >= 2.",
           call. = FALSE)
  }

  if (is.null(measurement.noise)) {

    mns <- 0

  } else if (is.numeric(measurement.noise)) {

    if (length(measurement.noise) != 1) {
        stop("`measurement.noise` must be a single value or a spectral object.",
             call. = FALSE)
    }

    # scale measurement noise variance by 2 * Nyquist frequency to get PSD level
    f.n <- max(spectra$mean$freq)
    mns <- c(measurement.noise) / (2 * f.n)

  } else {

    check.if.spectrum(measurement.noise)

    if (has.common.freq(measurement.noise, spectra$mean)) {
      measurement.noise <- InterpolateSpectrum(measurement.noise, spectra$mean)
    } else {
      stop("No sufficient frequency axis overlap between proxy data ",
           "and measurement noise spectrum.", call. = FALSE)
    }

    mns <- measurement.noise$spec

  }

  if (is.null(diffusion)) {

    dtf.corr <- 1

  } else {

    check.if.spectrum(diffusion)

    if (has.common.freq(diffusion, spectra$mean)) {
      diffusion <- InterpolateSpectrum(diffusion, spectra$mean)
    } else {
      stop("No sufficient frequency axis overlap between proxy data ",
           "and diffusion transfer function.", call. = FALSE)
    }

    dtf.corr <- 1 / diffusion$spec

  }

  if (is.null(time.uncertainty)) {

    ttf.corr <- 1

  } else {

    check.if.spectrum(time.uncertainty)

    if (has.common.freq(time.uncertainty, spectra$mean)) {
      time.uncertainty <- InterpolateSpectrum(time.uncertainty, spectra$mean)
    } else {
      stop("No sufficient frequency axis overlap between proxy data ",
           "and time uncertainty transfer function.", call. = FALSE)
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
  noise$spec  <- eta * dtf.corr * (mean - ttf.corr * stack - mns / eta)

  snr$spec <- signal$spec / noise$spec
  
  # organize output

  class(signal) <- class(noise) <- class(snr) <- "spec"

  setAttr <- function(x) {
    attr(x, "array.par") <- c(nc = N,
                              nt = attr(spectra, "array.par")[["nt"]],
                              res = attr(spectra, "array.par")[["res"]])
    return(x)
  }
  
  list(
    signal  = signal,
    noise   = noise,
    snr     = snr
  ) %>%
    setAttr()

}
