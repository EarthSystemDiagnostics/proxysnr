#' Diffusion transfer function
#'
#' This function implements an empirical Monte Carlo approach to estimate the
#' spectral transfer function for the effect of firn diffusion on the spatial
#' average of firn/ice-core stable isotope records.
#'
#' The approach is described in detail in Münch and Laepple (2018). In brief,
#' \code{nc} Gaussian white noise time series are created and diffused and the
#' average of these time series is calculated. The process is repeated
#' \code{ns} times. For each of the \code{ns} realisations, spectra of the
#' average diffused and undiffused records are calculated; subsequently, the
#' \code{ns} spectra are averaged, and the ratio of the average diffused to
#' the average undiffused spectrum yields the spectral transfer function.
#'
#' Diffusion is modelled as the convolution of the undiffused record with a
#' Gaussian with standard deviation given by the diffusion length
#' \code{sigma}. The spectral estimates are calculated using Thomson’s
#' multitaper method with three windows with linear detrending before
#' analysis.
#'
#' @param nt the length of the modelled isotope records (i.e. the number of
#'   data points in each record).
#' @param nc the number of cores in the modelled core array.
#' @param ns the number of Monte Carlo simulations for estimating the transfer
#'   function.
#' @param sigma numeric vector of length \code{nt} or numeric array of
#'   dimension \code{nt * nc} providing diffusion length values. The \code{nt}
#'   diffusion length values are assumed to correspond to the respective
#'   \code{nt} isotope values. If only a numeric vector is provided, it is
#'   assumed that these diffusion lengths are valid for all \code{nc} cores. If
#'   an array is provided, each column provides the diffusion lengths for the
#'   respective core. Note that the units of the diffusion length must match the
#'   units of \code{res}.
#' @param res the sampling (e.g., temporal) resolution of the isotope data;
#'   determines the frequency axis of the transfer function.
#' @param window length-2 vector giving a start and an end time (within
#'   \code{1:nt}) offering the possibility to only use a subset of the total
#'   length of the simulated records for the transfer function analysis, while
#'   the default of \code{NULL} means to use the records' entire lengths.
#' @param coherent if \code{TRUE}, \code{nc} identical white noise time series
#'   are assumed to estimate the transfer function; else (the default) \code{nc}
#'   independent noise series.
#' @param df.log width of the Gaussian kernel in logarithmic frequency units to
#'   smooth the spectral estimates; \code{NULL} (the default) suppresses
#'   smoothing. In general, smoothing should not be necessary when a
#'   sufficiently large number of Monte Carlo simulations (parameter \code{ns})
#'   is used; nevertheless, it can be switched on here when the transfer
#'   function still appears too noisy.
#' @param verbose.output logical controlling the size of the return object; per
#'   default, only the transfer function spectrum is returned, else also the
#'   spectra whose ratio determines the transfer function (see Details).
#' @param ... additional parameters which are passed to the spectral estimation
#'   function \code{\link{SpecMTM}}.
#'
#' @return either a spectral object (\code{?spec.object}) of the transfer function if
#'   \code{verbose.output = FALSE} (default), or a list of the spectral objects
#'   \code{signal}, \code{diffused} and \code{ratio}, providing the averages
#'   over the \code{ns} simulations of:
#'   \describe{
#'   \item{\code{signal}:}{the undiffused spectrum;}
#'   \item{\code{diffused}:}{the diffused spectrum;}
#'   \item{\code{ratio}}{their ratio (diffused/undiffused), i.e. the transfer
#'     function.}
#' }
#'
#' @seealso \code{\link{spec.object}} for the definition of a \code{proxysnr}
#'   spectral object.
#' @author Thomas Münch
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#'   2018.
#'
#' @export
#'
CalculateDiffusionTF <- function(nt, nc, ns, sigma, res = 1, window = NULL,
                                 coherent = FALSE, df.log = NULL,
                                 verbose.output = FALSE, ...) {

  # convert sigma vector into array if necessary and check for dimensions
  if (is.null(dim(sigma))) {
    
    if (length(sigma) != nt) {
      stop("Invalid length of supplied vector of diffusion lengths.")
    }
    
    sigma <- array(sigma, dim = c(nt, nc))

  } else {

    if (prod(dim(sigma)) != (nt * nc)) {
      stop("Invalid dimensions of supplied array of diffusion lengths.")
    }

  }

  if ((nw <- length(window)) > 0) {

    if (nw != 2) stop("`window` must have length two.", call. = FALSE)
    if (window[2] <= window[1])
      stop("`window[2]` must be > `window[1]`.", call. = FALSE)
    if (window[1] < 1) stop("`window[1]` must be >= 1.", call. = FALSE)
    if (window[2] > nt)
      stop("`window[2]` is > total number of time points.", call. = FALSE)

    k <- window[1] : window[2]

  } else {

    k <- seq_len(nt)

  }

  if ((nsmooth <- length(df.log)) > 1)
    stop("`df.log` must be of length 1 or `NULL`.")

  apply.smoothing <- nsmooth == 1

  # function for forward diffusion
  .diffuse <- function(rec, sigma, res = 1){

    n <- length(rec)

    # scale diffusion length according to resolution of record
    sigma <- sigma / res

    # pad end of record with mean of record to avoid NA's
    # at the end of diffused record
    rec <- c(rec, rep(mean(rec, na.rm = TRUE), 10 * max(sigma)))

    # vector to store diffused data
    rec.diffused <- rep(NA, n)

    # loop over record
    for (i in 1 : n) {

      # diffusion length for current depth
      sig <- sigma[i]

      # set range of convolution integral (= 2*imax + 1) to ~ 10*sig
      imax <- ceiling(5 * sig)
      ran <- (i - imax) : (i + imax)

      # skip part of range in convolution integral which extends above surface
      ran <- ran[ran > 0]
      # relative range for convolution kernel
      rel.ran <- i - ran

      # convolution kernel
      kernel <- exp(-(rel.ran)^2 / (2 * sig^2))
      kernel <- kernel / sum(kernel)

      # diffuse data at current depth bin
      diff.value <- sum(rec[ran] * kernel)

      rec.diffused[i] <- diff.value
    }

    return(rec.diffused)

  }

  # create 'nc' white noise time series, diffuse them, and calculate their
  # average; repeat 'ns' times
  
  stack.signal <- array(dim = c(nt, ns))
  stack.diff   <- array(dim = c(nt, ns))

  for (i in 1 : ns) {

    if (coherent) {
      X <- array(stats::rnorm(nt), dim = c(nt, nc))
    } else {
      X <- array(stats::rnorm(nt * nc), dim = c(nt, nc))
    }
    
    Xdiff <- array(dim = c(nt, nc))

    for (j in 1 : nc) {

      Xdiff[, j] <- .diffuse(X[, j], sigma = sigma[, j], res = res)
    }

    stack.signal[, i] <- rowMeans(X)
    stack.diff[, i]   <- rowMeans(Xdiff)

  }


  # calculate spectra of the 'ns' average diffused and undiffused noise series

  signal.spec <- lapply(seq_len(ncol(stack.signal)), function(i) {
    SpecMTM(stats::ts(stack.signal[k, i], deltat = res), ...)})

  diff.spec <- lapply(seq_len(ncol(stack.diff)), function(i) {
    SpecMTM(stats::ts(stack.diff[k, i], deltat = res), ...)})


  # calculate the average over all 'ns' simulations

  signal.spec.mean <- MeanSpectrum(signal.spec)
  diff.spec.mean   <- MeanSpectrum(diff.spec)


  # return average undiffused and diffused spectra and the corresponding
  # ratio (transfer function) as a list

  version <- sprintf("Creation date: %s.", Sys.time())
  nsim <- sprintf("Number of simulations used: N = %s.",
                  formatC(ns, big.mark = ",", format = "d"))
  smoothing <- paste("Log-smooth applied:",
    if (apply.smoothing) sprintf("Yes (df.log = %1.2f).", df.log) else "No.")

  setAttr <- function(x) {
    attr(x, "version") <- version
    attr(x, "N.sim") <- nsim
    attr(x, "log-smooth") <- smoothing
    return(x)
  }
  
  res <- list(
  signal   = signal.spec.mean,
  diffused = diff.spec.mean,
  ratio    = list(freq = signal.spec.mean$freq,
                  spec = diff.spec.mean$spec / signal.spec.mean$spec)
  ) %>%
    {if (apply.smoothing) {lapply(., LogSmooth, df.log = df.log)} else {.}} %>%
    lapply(setAttr)

  class(res$ratio) <- "spec"

  if (verbose.output) return(res) else return(res$ratio)

}
