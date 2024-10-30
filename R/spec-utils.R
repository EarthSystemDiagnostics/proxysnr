# Functions in this file are copied from R package paleospec:
# https://github.com/EarthSystemDiagnostics/paleospec
# 
# Original authors and copyright information follow:
# 
# * LogSmooth, SpecMTM, LPlot, LLines:
#   Written by Thomas Laepple. MIT, Copyright (C) 2019 Earth System Diagnostics
#   group of the Alfred Wegener Institute
# * MeanSpectrum:
#   Written by Thomas Laepple, Thomas Münch. MIT, Copyright (C) 2019 Earth
#   System Diagnostics group of the Alfred Wegener Institute
# * SimPLS:
#   Written by Andrew M Dolman, Torben Kunz. MIT, Copyright (C) 2019 Earth
#   System Diagnostics group of the Alfred Wegener Institute

#' Spectral smoothing
#'
#' Smooth spectrum using Gaussian kernel with constant width in logarithmic
#' frequency space.
#'
#' @param x a spectral object including the degrees of freedom of the spectrum.
#' @param df.log smoothing width in logarithmic frequency units.
#'
#' @return a spectral object with the smoothed spectrum including the average
#'   degrees of freedom.
#'
#' @author Thomas Laepple
#' @noRd
#'
LogSmooth <- function(x, df.log=0.05) {

  # Gaussian kernel
  weights <- function(x, sigma) {
    1 / sqrt(2 * pi * sigma^2) * exp(-x^2 / (2 * sigma^2))
  }

  # frequency dependent smoothing weights in log space
  fweights <- function(ftarget, f, df.log) {
    
    sigma   <- ftarget * (exp(df.log) - exp(-df.log))
    weights <- weights(f - ftarget, sigma)
    
    return(weights)
  }

  # logarithmic smoothing with cut weights, if necessary
  smoothlog.cutEnd <- function(x, f, df.log, dof = 1) {

    x.smooth <- vector()
    dof.smooth <- vector()
    for (i in 1 : length(f)) {

      # get averaging weights
      w <- fweights(f[i], f, df.log)
      # cut the weights down
      DistanceSlowEnd <- i - 1
      DistanceFastEnd <- length(f) - i

      if ((i + DistanceSlowEnd + 1) <= length(f))
        w[(i + DistanceSlowEnd + 1) : length(f)] <- 0
      if ((i - DistanceFastEnd - 1) >= 1)
        w[1 : (i - DistanceFastEnd - 1)] <- 0

      # normalize to 1
      w <- w / f
      w <- w / sum(w)

      # apply the smoothing
      x.smooth[i]   <- sum(w * x)
      dof.smooth[i] <- sum(w * dof) / sum(w^2)
    }
    
    return(list(spec = x.smooth, dof = dof.smooth))
    
  }

  temp <- smoothlog.cutEnd(x$spec, x$freq, df.log, dof = x$dof)

  result <- list()
  
  result$freq <- x$freq
  result$spec <- temp$spec
  result$dof <- temp$dof

  class(result) <- "spec"
  
  return(result)

}

#' Multitaper spectral estimate
#' 
#' This is a wrapper for the \code{\link[multitaper]{spec.mtm}} function from
#' package multitaper, which sets convenient defaults for the spectral estimate
#' and makes the degrees of freedom accessible as a direct list element of the
#' output.
#'
#' @inheritParams multitaper::spec.mtm
#' @param detrend Shall the time series data be linearly detrended before
#'   computing the spectrum? Defaults to \code{TRUE}.
#' @param bPad Shall the time series data be padded before computing the
#'   spectrum? Defaults to \code{FALSE}.
#' @param ... further parameters passed to \code{\link[multitaper]{spec.mtm}}.
#' @return a spectral object including the degrees of freedom of the spectral
#'   estimate (copied from the \code{spec.mtm} default output).
#'
#' @author Thomas Laepple
#' @seealso \code{\link[multitaper]{spec.mtm}}
#'
#' @references
#'
#' Thomson, D.J (1982) Spectrum estimation and harmonic analysis. _Proceedings
#' of the IEEE_ Volume *70*, Number 9, pp. 1055-1096.
#'
#' Percival, D.B. and Walden, A.T. (1993) _Spectral analysis for physical
#' applications_ Cambridge University Press.
#'
SpecMTM <- function(timeSeries, k = 3, nw = 2, nFFT = "default",
                    centre = c("Slepian"), dpssIN = NULL,
                    returnZeroFreq = FALSE, Ftest = FALSE, jackknife = FALSE,
                    jkCIProb = 0.95, maxAdaptiveIterations = 100,
                    plot = FALSE, na.action = stats::na.fail,
                    returnInternals = FALSE, detrend = TRUE,
                    bPad = FALSE, ...) {

  if (sum(is.na(timeSeries)) > 0)
    stop("missing data")
  if (!bPad)
    nFFT = length(timeSeries)
  if (detrend)
    timeSeries[] <- stats::lm(timeSeries ~ seq(timeSeries))$residuals

  result <- multitaper::spec.mtm(timeSeries = timeSeries,
                                 k = k, nw = nw, nFFT = nFFT,
                                 centre = centre,
                                 dpssIN = dpssIN,
                                 returnZeroFreq = returnZeroFreq,
                                 Ftest = Ftest,
                                 jackknife = jackknife,
                                 jkCIProb = jkCIProb,
                                 maxAdaptiveIterations = maxAdaptiveIterations,
                                 plot = plot,
                                 na.action = na.action,
                                 returnInternals = returnInternals,
                                 ...)
  result$dof <- result$mtm$dofs
  
  return(result)
}

#' Mean spectrum
#'
#' Calculate the mean spectrum of all supplied individual spectra.
#'
#' The mean spectrum is calculated as a simple average across all power spectral
#' densities from the individual spectra. Degrees of freedom (\code{dof}) of the
#' mean spectrum are obtained from the sum of the individual dof's at each
#' frequency.
#'
#' @param speclist a list of spectral objects.
#'
#' @return a spectral object with the averaged spectrum.
#'
#' @author Thomas Laepple, Thomas Münch
#' @noRd
#'
MeanSpectrum <- function(speclist) {

  # check for equal lengths of supplied spectra
  if (stats::var(lengths(lapply(speclist, "[[", "freq"))) > 0)
    stop("MeanSpectrum: Spectra are of different lengths.", call. = FALSE)
  
  mean <- list()
  mean$freq <- speclist[[1]]$freq
  mean$spec <- rowMeans(sapply(speclist, function(x) {x$spec}))
  mean$dof  <- rowSums(sapply(speclist, function(x) {x$dof}))

  class(mean) <- "spec"

  return(mean)

}

#' Log-log spectral plot
#'
#' This function plots a spectrum on a double-logarithmic scale and optionally
#' adds a transparent confidence interval.
#'
#' @param x a spectral object.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#'   interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} the x-axis is displayed in units of period
#'   (inverse frequency), increasing to the left. Defaults to \code{FALSE}.
#' @param bNoPlot if \code{TRUE} only produce the plot frame (\code{type = "n"}
#'   behaviour of function \code{\link{plot}}). Defaults to \code{FALSE}.
#' @param axes if \code{FALSE} the plotting of the x and y axes is
#'   suppressed. Defaults to \code{TRUE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#'   interval. Defaults to \code{0.2}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side.
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param xlab character string for labelling the x-axis.
#' @param ylab character string for labelling the y-axis.
#' @param xlim range of x-axis values; if \code{NULL} (the default) it is
#'   calculated internally and automatically reversed for \code{bPeriod = TRUE}.
#' @param ylim range of y-axis values; if \code{NULL} (the default) it is
#'   calculated internally.
#' @param ... further graphical parameters passed to \code{plot}.
#'
#' @author Thomas Laepple
#' @noRd
#'
LPlot <- function(x, conf = TRUE, bPeriod = FALSE, bNoPlot = FALSE, axes = TRUE,
                  col = "black", alpha = 0.2, removeFirst = 0, removeLast = 0,
                  xlab = "f", ylab = "PSD", xlim = NULL, ylim = NULL, ...) {

  if (bPeriod) {
    x$freq <- 1 / x$freq
    if (is.null(xlim)) xlim <- rev(range(x$freq))
    if (xlab == "f") xlab <- "period"
  }

  index <- (removeFirst + 1) : (length(x$freq) - removeLast)
  
  x$freq  <- x$freq[index]
  x$spec  <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]

  graphics::plot(x$freq, x$spec, type = "n", log = "xy",
                 xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                 axes = axes, ...)

  lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
  if (conf & lim & !bNoPlot) {
    Polyplot(x = x$freq, y1 = x$lim.1, y2 = x$lim.2,
             col = col, alpha = alpha)
  }

  if (!bNoPlot) graphics::lines(x$freq, x$spec, col = col, ...)
  
}

#' Add spectrum to existing log-log spectral plot
#'
#' This function adds a spectrum to an existing double-logarithmic plot and
#' optionally adds a transparent confidence interval.
#' 
#' @param x a spectral object.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#'   interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} treat the x-axis values in units of period
#'   (inverse frequency). Defaults to \code{FALSE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#'   interval. Defaults to \code{0.2}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side. 
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param ... further graphical parameters passed to \code{lines}.
#'
#' @author Thomas Laepple
#' @noRd
#'
LLines<-function(x, conf = TRUE, bPeriod = FALSE, col = "black", alpha = 0.2,
                 removeFirst = 0, removeLast = 0, ...) {

  if (bPeriod) x$freq <- 1 / x$freq
  
  index <- (removeFirst + 1) : (length(x$freq) - removeLast)

  x$freq  <- x$freq[index]
  x$spec  <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]

  lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
  if (conf & lim) {
    Polyplot(x = x$freq, y1 = x$lim.1, y2 = x$lim.2,
             col = col, alpha = alpha)
  }
  
  graphics::lines(x$freq, x$spec, col = col, ...)
  
}

#' Simulate a random time series with a power-law spectrum
#'
#' This function creates a power-law series. It has the problem that it
#' effectively produces (fractional) Brownian bridges, that is, the end is close
#' to the beginning (cyclic signal), rather than true fBm or fGn series.
#'
#' @param N length of time series to be generated.
#' @param beta slope of the power law: `beta = 1` produces time series with a
#'   slope of -1 when plotted on log-log power ~ frequency axes.
#' @param alpha a constant. If `alpha > 0` this is the power-law parameter in
#'   `alpha * f^(-beta)`. If `alpha < 0`, the variance of the returned
#'   time series is scaled so that its expected value is `abs(alpha)` and its
#'   expected PSD is proportional to `f^(-beta)`.
#' @return a vector of length `N`.
#'
#' @author Torben Kunz, Andrew Dolman
#' @noRd
#'
SimPLS <- function(N, beta = 0, alpha = -1) {

  frequency.axis <- function(N) {
    fax <- 0 : (N - 1) / N
    fax[fax > 0.5] <- fax[fax > 0.5] - 1
    fax
  }

  # Pad the length of the timeseries so that it is highly composite - this speeds
  # up the FFT operations.
  N2 <- (3^ceiling(log(N, base = 3)))

  x2 <- stats::rnorm(N2)
  xfft <- stats::fft(x2)

  fax <- frequency.axis(N2)

  P <- c(0, abs(alpha) * abs(fax[2 : N2])^(-beta))

  if (alpha < 0) P <- P * abs(alpha) * (N2 - 1) / sum(P[1 : N2])

  xfft <- xfft * sqrt(P)
  x2 <- stats::fft(xfft, inverse = TRUE) / N2
  x2 <- Re(x2)

  # trim the timeseries back to the requested length
  x <- x2[1 : N]

  # scale the variance of timeseries at requested length
  if (alpha < 0) {
    sdx2 <- stats::sd(x2)
    x <- sdx2 * x / stats::sd(x)
  }

  return(x)
}
