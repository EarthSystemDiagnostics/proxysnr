#' Calculate signal and noise spectra
#' 
#' Calculate the raw signal and noise spectra, and the corresponding
#' signal-to-noise ratio, as estimated from a core array of \code{n} proxy
#' records, and correct these, where applicable, for the effects of time
#' uncertainty and water vapour diffusion (relevant for firn and ice cores).
#'
#' This function is an implementation of Eq. (4) in Münch and Laepple
#' (2018). While the diffusion correction is relevant only for diffusing
#' proxies such as stable isotopes from firn and ice cores, this function can
#' be applied to a large set of proxy data since only one of the two
#' correction functions, or none, need to be supplied; thus, e.g., it can also
#' be used for non-diffusing proxy data where only time uncertainty is
#' relevant, or for estimating raw signal and noise spectra by supplying no
#' correction functions at all.
#'
#' @param spectra a list of the raw spectral estimates from a proxy core
#'   array. Expected is the output from \code{\link{ObtainArraySpectra}}, but
#'   sufficient is a named list of two components giving the \code{mean} and
#'   \code{stack} spectra.
#' @param neff the effective number of records (e.g. to account for an expected
#'   spatial correlation of the local noise). Per default set to element
#'   \code{N} in \code{spectra}, otherwise supply it explicitly here.
#' @param diffusion numeric vector of diffusion correction values (inverse
#'   transfer function); must be of the same length as the spectral estimates in
#'   \code{spectra}. The default `NULL` is to apply no correction.
#' @param time.uncertainty numeric vector of time uncertainty correction
#'   values (inverse transfer function); must be of the same length as the
#'   spectral estimates in \code{spectra}. The default `NULL` is to apply no
#'   correction.
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
#' @seealso \code{\link{ObtainArraySpectra}}
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

  if (!all(hasName(spectra, c("mean", "stack"))))
    stop("`spectra` must have elements `mean` and `stack`.", call. = FALSE)

  is.spectrum(spectra$mean)
  is.spectrum(spectra$stack)

  if (length(spectra$mean$freq) != length(spectra$stack$freq)) {
    stop("`mean` and `stack` must have the same number of spectral estimates.",
         call. = FALSE)
  }

  if (is.null(neff)) {
    stop("Supply (effective) number of records.")
  }

  if (is.null(diffusion)) {
    diffusion <- 1
  } else {
    if (length(diffusion) != length(spectra$mean$freq)) {
      stop("Length of diffusion correction ",
           "does not match length of spectral estimates.")
    }
  }

  if (is.null(time.uncertainty)) {
    time.uncertainty <- 1
  } else {
    if (length(time.uncertainty) != length(spectra$mean$freq)) {
      stop("Length of time uncertainty correction ",
           "does not match length of spectral estimates.")
    }
  }


  # calculate the signal and noise spectra
  
  d.corr <- diffusion
  t.corr <- time.uncertainty

  N <- neff
  n <- N / (N - t.corr)

  corr.fac <- n * d.corr

  mean  <- spectra$mean
  stack <- spectra$stack
  
  signal <- list()
  noise  <- list()
  snr    <- list()

  signal$freq <- mean$freq
  noise$freq  <- mean$freq
  snr$freq    <- mean$freq

  signal$spec <- t.corr * corr.fac * (stack$spec - mean$spec / N)
  noise$spec  <- corr.fac * (mean$spec - t.corr * stack$spec)

  snr$spec <- signal$spec / noise$spec

  
  # Organize output

  class(signal) <- "spec"
  class(noise)  <- "spec"
  class(snr)    <- "spec"
  
  res <- list(
    signal  = signal,
    noise   = noise,
    snr     = snr)

  return(res)

}
