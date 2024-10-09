#' Simulate a proxy core array
#'
#' Simulate an array of proxy records with each record being composed of a
#' signal that is common to all the records in the array and an independent
#' noise component. Both components are simulated such that their spectra follow
#' a power-law of the form `alpha * f^(-beta)` where `f` denotes frequency.
#'
#' @param signal.par a named list with elements `alpha` and `beta` supplying the
#'   power-law parameters to simulate the common signal.
#' @param noise.par as \code{signal.par} the power-law parameters to simulate
#'   the independent noise components.
#' @param nc the number of proxy records (`"cores"`) in the array.
#' @param nt the length of each record in the array.
#' @return a list with `nc` elements where each is one simulated proxy record of
#'   length `nt`.
#'
#' @author Thomas Münch
#'
simCoreArray <- function(signal.par, noise.par, nc, nt) {

  # simulate powerlaw signal
  signal <- SimPLS(N = nt, beta = signal.par$beta, alpha = signal.par$alpha)

  # simulate nc proxy records (= signal + independent noise)
  cores <- lapply(seq(nc), function(i) {

    signal + SimPLS(N = nt, beta = noise.par$beta, alpha = noise.par$alpha)

  })

  return(cores)

}

#' Simulate a core array and estimate the signal and noise spectra
#'
#' This is a wrapper function to first simulate a proxy core array
#' (\code{\link{simCoreArray}}), estimate the corresponding array spectra
#' (\code{\link{ObtainArraySpectra}}), and finally separate for these simulated
#' data the common signal from the independent noise
#' (\code{\link{SeparateSignalFromNoise}}).
#'
#' @param res the sampling (e.g., temporal) resolution of the simulated proxy
#'   records; determines the frequency axis of the spectral estimates.
#' @param df.log width of the Gaussian kernel in logarithmic frequency units to
#'   smooth the spectral estimates; \code{NULL} (the default) suppresses
#'   smoothing.
#' @inheritParams simCoreArray
#' @return a named list of the three spectral objects `signal`, `noise` and
#'   `snr`, giving the spectra of the estimated common signal, the independent
#'   noise, and the corresponding signal-to-noise ratio (SNR), respectively.
#'
#' @author Thomas Münch
#'
simSignalAndNoise <- function(signal.par, noise.par, nc, nt, res,
                              df.log = NULL) {

  simCoreArray(signal.par, noise.par, nc, nt) %>%
    ObtainArraySpectra(res = res, df.log = df.log) %>%
    SeparateSignalFromNoise()

}

#' Monte Carlo signal and noise simulations 
#'
#' Replicate \code{\link{simSignalAndNoise}} a given number of times.
#'
#' @param nmc integer; the number of replications.
#' @inheritParams simSignalAndNoise
#' @return a list of length \code{nmc} where each list element is one signal and
#'   noise estimation realization, i.e. the output of
#'   \code{\link{simSignalAndNoise}}.
#' @author Thomas Münch
#' @seealso \code{\link{simSignalAndNoise}}, \code{\link{simCoreArray}}
#'
runSurrogates <- function(signal.par, noise.par, nc, nt, res, nmc = 10,
                          df.log = NULL) {

  replicate(nmc,
            simSignalAndNoise(signal.par = signal.par,
                              noise.par = noise.par,
                              nc = nc, nt = nt, res = res,
                              df.log = df.log),
            simplify = FALSE)

}
