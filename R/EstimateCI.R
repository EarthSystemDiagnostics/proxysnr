#' Confidence interval estimation
#'
#' This function estimates based on a parametric bootstrapping procedure (see
#' Details) confidence intervals for a set of signal, noise and SNR spectra
#' obtained from an actual proxy record array.
#'
#' The parametric bootstrapping procedure for the confidence level estimation is
#' implemented as follows. A power-law fit of the form \code{alpha * f^(-beta)}
#' is applied to the actual signal and noise spectra, and the resulting
#' power-law coefficients are used to generate surrogate signal and noise series
#' in a simulated array that mimics the actual proxy record array. This
#' simulated array is replicated \code{nmc} times and lower and upper quantiles
#' are calculated across the realizations of the signal, noise and SNR
#' surrogates. Subsequently, the quantiles are scaled to the respective mean
#' estimates and applied multiplicatively to the actual proxy estimates of
#' signal, noise, and SNR to yield the confidence intervals.
#'
#' @param spectra a list with the spectral objects \code{signal}, \code{noise},
#'   and \code{snr} for an investigated proxy record array, obtained from
#'   \code{\link{SeparateSignalFromNoise}}.
#' @param f.start lower end of the frequency range on which the power-law fit is
#'   made on the proxy data (see Details); the default \code{NULL} uses the
#'   lowest frequency of the proxy \code{spectra}.
#' @param f.end as \code{f.start} for the upper end; the default \code{NULL}
#'   uses the uppermost frequency of the proxy \code{spectra}.
#' @param nmc integer; the number of replications for the confidence interval
#'   estimation.
#' @param probs length-2 numeric vector of probabilities with values in [0,1]
#'   defining the confidence interval; defaults to the 10-90 \% interval.
#' @param df.log width of the Gaussian kernel in logarithmic frequency units to
#'   smooth the spectral estimates of the simulated data; some smoothing is
#'   usually necessary to avoid physically implausible negative power
#'   occasionally occuring for some frequencies upon estimating the common
#'   signal spectrum. It is suggested to use the same amount of smoothing as for
#'   the actual proxy data, while setting \code{NULL} (the default) suppresses
#'   smoothing.
#' @param ci.df.log width of the Gaussian smoothing kernel to smooth the
#'   estimated confidence intervals, merely for visual purposes; \code{NULL}
#'   (the default) suppresses smoothing.
#' @return the input \code{spectra} object, amended by the confidence intervals
#'   for the signal, noise and SNR spectra (element \code{lim.1} gives the upper
#'   confidence level, element \code{lim.2} the lower level, respectively).
#'
#' @author Thomas Münch
#' @seealso \code{\link{SeparateSignalFromNoise}},
#'   \code{\link{ObtainArraySpectra}}
#'
#' @export
#'
EstimateCI <- function(spectra, f.start = NULL, f.end = NULL, nmc = 10,
                       probs = c(0.1, 0.9), df.log = NULL, ci.df.log = NULL) {

  # ----------------------------------------------------------------------------
  # input error checking

  if (!is.list(spectra)) stop("`spectra` must be a list.", call. = FALSE)

  if (!all(utils::hasName(spectra, c("signal", "noise", "snr"))))
    stop("`spectra` must have elements `signal`, `noise`, and `snr`.",
         call. = FALSE)

  check.if.spectrum(spectra$signal)
  check.if.spectrum(spectra$noise)
  check.if.spectrum(spectra$snr)

  if (length(spectra$signal$freq) != length(spectra$noise$freq)) {
    stop("`signal` and `noise` must have the same number of spectral estimates.",
         call. = FALSE)
  }
  if (length(spectra$signal$freq) != length(spectra$snr$freq)) {
    stop("`signal` and `snr` must have the same number of spectral estimates.",
         call. = FALSE)
  }

  if (!all(spectra$signal$freq == spectra$noise$freq)) {
    stop("Frequency axes of `signal` and `noise` do not match.", call. = FALSE)
  }
  if (!all(spectra$signal$freq == spectra$snr$freq)) {
    stop("Frequency axes of `signal` and `snr` do not match.", call. = FALSE)
  }

  has.array.attribute(spectra)

  # ----------------------------------------------------------------------------
  # helper functions

  # set limits on spectral object multiplicatively from estimated CI
  .setlimits <- function(x, lower, upper) {

    x$lim.1 <- upper * x$spec
    x$lim.2 <- lower * x$spec

    return(x)

  }

  # log-smooth estimated CI
  .smooth <- function(ci, df.log = 0.05) {

  if (is.null(df.log)) return(ci)

  ci$lower <- LogSmooth(list(freq = ci$freq, spec = ci$lower),
                        df.log = df.log)$spec
  ci$upper <- LogSmooth(list(freq = ci$freq, spec = ci$upper),
                        df.log = df.log)$spec

  return(ci)

  }

  # interpolate estimated CI
  .interpolate.ci <- function(ci, target) {

    suppressWarnings(
      ci$lower <- InterpolateSpectrum(list(freq = ci$freq, spec = ci$lower),
                                      target)$spec
    )
    suppressWarnings(
      ci$upper <- InterpolateSpectrum(list(freq = ci$freq, spec = ci$upper),
                                      target)$spec
    )

    return(ci)

  }

  # check if simulated frequency axis overlaps with input frequency axis and
  # interpolate if needed
  .check.frequency.axes <- function(sim, target) {

    if (!has.common.freq(sim$signal, target)) {

      warning("Input and simulated frequency axes differ;\n",
              "confidence intervals are interpolated possibly ",
              "producing NA values.",
              call. = FALSE)

      lapply(sim, .interpolate.ci, target = target)

    } else {

      sim

    }

  }

  # ----------------------------------------------------------------------------
  # run CI estimation

  # run surrogate signal and noise simulation, extract confidence intervals,
  # optionally smooth them, and interpolate to input frequency axis if needed

  ci <- runSimulation(spectra = spectra, f.start = f.start, f.end = f.end,
                      nmc = nmc, df.log = df.log) %>%
    extractQuantiles(probs = probs) %>%
    lapply(.smooth, df.log = ci.df.log) %>%
    .check.frequency.axes(spectra$signal)

  # supply quantiles to input data and return

  spectra$signal <- .setlimits(spectra$signal, ci$signal$lower, ci$signal$upper)
  spectra$noise  <- .setlimits(spectra$noise, ci$noise$lower, ci$noise$upper)
  spectra$snr    <- .setlimits(spectra$snr, ci$snr$lower, ci$snr$upper)

  return(spectra)

}

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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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

#' Run simulation for the CI estimation
#'
#' Get the required parameters from the input spectra (number of observations,
#' observation resolution, number of proxy records, and signal and noise
#' powerlaw parameters) and call the Monte Carlo signal and noise simulation
#' \code{\link{runSurrogates}}.
#'
#' @inheritParams EstimateCI
#' @return the output of \code{runSurrogates}: a list of length \code{nmc}
#'   where each list element is one signal and noise estimation realization.
#' @author Thomas Münch
#' @seealso \code{\link{runSurrogates}}, \code{\link{EstimateCI}}
#' @keywords internal
#'
runSimulation <- function(spectra, f.start = NULL, f.end = NULL,
                          nmc = 10, df.log = NULL) {

  # estimate powerlaw coefficients
  signal.par <- fit.powerlaw(spectra$signal, f.start, f.end)
  noise.par  <- fit.powerlaw(spectra$noise, f.start, f.end)

  # get number of cores underlying results in `spectra`
  nc <- attr(spectra, "array.par")[["nc"]]
  # get number of time steps of original proxy records
  nt <- attr(spectra, "array.par")[["nt"]]
  # get time resolution of original proxy records
  res <- attr(spectra, "array.par")[["res"]]

  runSurrogates(signal.par = signal.par, noise.par = noise.par,
                nc = nc, nt = nt, res = res, nmc = nmc,
                df.log = df.log)

}

#' Quantiles of simulated signal and noise spectrum realizations
#'
#' Extract lower and upper quantiles across the signal, noise and SNR
#' realizations from a run of \code{\link{runSurrogates}}, scaled to the
#' respective mean estimates of the realiztaions.
#'
#' @param surrogates the output from a run of \code{\link{runSurrogates}}.
#' @param probs length-2 numeric vector of probabilities with values in [0,1];
#'   default setting extracts the 10 and 90 \% quantiles.
#' @param minval minimum value assigned to negative estimates of the scaled
#'   lower and upper quantiles to avoid unphysical negative confidence
#'   intervals, thereby preventing issues on logarithmic plots.
#' @return a named list of the three elements `signal`, `noise`, and `snr`, with
#'   each list element in turn being a named list with the `lower` and `upper`
#'   quantiles (according to the setting of \code{probs}) across the simulated
#'   signal, noise and signal-to-noise ratio (SNR) spectra.
#'
#' @author Thomas Münch
#' @seealso \code{\link{runSurrogates}}
#' @keywords internal
#'
extractQuantiles <- function(surrogates, probs = c(0.1, 0.9), minval = 1e-15) {

  if (length(probs) != 2)
    stop("`probs` must be a length-2 vector.", call. = FALSE)

  if (diff(probs) <= 0)
    stop("`probs[2]` must be > `probs[1]`.", call. = FALSE)

  .extract.runs <- function(run, name) run[[name]][["spec"]]

  .extract.quantiles <- function(surrogates, name) {

    f <- surrogates[[1]][["signal"]][["freq"]]
    x <- sapply(surrogates, .extract.runs, name = name)

    m <- rowMeans(x)
    m[m <= 0] <- minval

    lower <- apply(x, 1, stats::quantile, probs = probs[1]) / m
    upper <- apply(x, 1, stats::quantile, probs = probs[2]) / m

    upper[upper <= 0] <- minval
    lower[lower <= 0] <- minval

    list(freq = f, lower = lower, upper = upper)

  }

  list(
    signal = .extract.quantiles(surrogates, "signal"),
    noise  = .extract.quantiles(surrogates, "noise"),
    snr    = .extract.quantiles(surrogates, "snr"))

}
