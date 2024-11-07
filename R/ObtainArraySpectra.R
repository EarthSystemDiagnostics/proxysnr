#' Spectral estimates of core array
#'
#' Calculate all relevant spectral estimates for a given array of \code{n}
#' proxy records. The spectral estimates can be smoothed in logarithmic space
#' and are calculated using Thomson’s multitaper method with three windows with
#' linear detrending before analysis.
#'
#' @param cores a list or a data frame of the proxy data from the core array. If
#'   a list is supplied, all elements must be numeric vectors of the same
#'   length.
#' @param res the sampling (e.g., temporal) resolution of the proxy data;
#'   determines the frequency axis of the spectral estimates.
#' @param neff the effective number of records (e.g. to account for an expected
#'   spatial correlation of the local noise). Per default, no spatial
#'   correlation is assumed and \code{neff} is set to the number of proxy
#'   records (the length of \code{cores}).
#' @param df.log width of the Gaussian kernel in logarithmic frequency units to
#'   smooth the spectral estimates; \code{NULL} (the default) suppresses
#'   smoothing.
#' @param ... additional parameters which are passed to the spectral estimation
#'   function \code{\link{SpecMTM}}.
#'
#' @return A list of the following components:
#'   \describe{
#'   \item{single:}{a list of \code{N} spectral objects (`?spec.object`) with
#'     the spectra of each individual proxy record;}
#'   \item{mean:}{spectral object of the mean spectrum across all individual
#'     spectra;}
#'   \item{stack:}{spectral object of the spectrum of the average proxy record
#'     in the time domain ("stacked record");}
#' }
#' with the attribute "array.par": a named vector with information on the proxy
#' record array: number of (effective) records ("nc" = \code{neff}), number of
#' observation points per record ("nt"), and sampling resolution ("res" =
#' \code{res}).
#'
#' @author Thomas Münch
#' @seealso \code{\link{PlotArraySpectra}}, \code{\link{SpecMTM}},
#'   `?spec.object` for the definition of a "proxysnr" spectral object.
#'
#' @export
#'
ObtainArraySpectra <- function(cores, res = 1, neff = length(cores),
                               df.log = NULL, ...) {

  # error checking
  if (!is.list(cores))
    stop("`cores` must be a list or a data frame.", call. = FALSE)

  if (length(cores) == 1)
    stop("Only one proxy record (`cores` is of length 1).", call. = FALSE)

  ll <- lengths(cores, use.names = FALSE)
  if (!is.data.frame(cores)) {

    # data vectors in list must be of the same length
    if (stats::sd(ll) > 0) {
      stop("Elements of `cores` must all have the same length.", call. = FALSE)
    }

  }

  nt <- ll[1]
  if (nt <= 8) {
    stop("Need at least 9 observations per record to calculate spectra.",
         call. = FALSE)
  }
  
  # estimate individual spectra
  single <- cores %>%
    lapply(stats::ts, deltat = res) %>%
    lapply(SpecMTM, ...)

  # estimate spectrum of stacked record
  stack <- simplify2array(cores) %>%
    rowMeans %>%
    stats::ts(deltat = res) %>%
    SpecMTM(...)

  # calculate mean spectrum across individual record's spectra
  mean <- MeanSpectrum(single)

  # log-smooth spectra
  if (!is.null(df.log)) {

    single <- lapply(single, LogSmooth, df.log = df.log)
    mean   <- LogSmooth(mean, df.log = df.log)
    stack  <- LogSmooth(stack, df.log = df.log)
  }

  # return results as a list with attribute supplying the array parameter

  setAttr <- function(x) {
    attr(x, "array.par") <- c(nc = neff, nt = nt, res = res)
    return(x)
  }

  list(
    single = single,
    mean   = mean,
    stack  = stack
  ) %>%
    setAttr()

}
