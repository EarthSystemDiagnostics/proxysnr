##' Frequency-dependent signal-to-noise ratio from integrated spectra
##'
##' This function calculates the signal-to-noise ratio (SNR) as a function of
##' frequency, interpreted as the temporal resolution of a proxy record. This
##' variant of the SNR is obtained from the ratio of a signal and a noise
##' spectrum, which are integrated before cumulatively across frequencies.
##'
##' The function is an implementation of Eq. (6) in Münch and Laepple
##' (2018). The integral in (6) is approximated by the cumulative sum of the
##' integration arguments from \code{f.int1} to \code{f.int2}, where
##' \code{f.int1 = f1} and \code{f.int2} consecutively increases from \code{f1}
##' to \code{f2}.
##' @param input a list of the spectral objects lists \code{signal} and
##' \code{noise}, usually to be obtained from a call to
##' \code{\link{SeparateSpectra}}
##' @param N integer; number of proxy records averaged. The default returns the
##' SNR assuming a single proxy record. For a different number, the SNR is
##' scaled by this number, assuming independent noise between the records.
##' @param f1 index of the the minimum frequency from which to integrate the
##' signal and noise spectra for calculating the SNR; per default the
##' lowest frequency of the spectral estimates is omitted
##' @param f2 index of the maximum frequency until which to integrate the signal
##' and noise spectra for calculating the SNR; defaults to use the
##' maximum frequency of the given spectral estimates
##' @inheritParams StackCorrelation
##' @return a spectral object list of the SNR.
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal- to centennial-scale isotope variations from Antarctic ice cores?
##' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
##' @export
IntegratedSNR <- function(input, N = 1, f1 = 2, f2 = "max",
                          freq.cut.lower = NULL, freq.cut.upper = NULL) {

  if (!is.list(input)) stop("'input' needs to be a list.", call. = FALSE)
  if (length(stats::na.omit(match(names(input), c("signal", "noise")))) != 2) {
    stop("'input' list needs to contain the elements 'signal' and 'noise'.",
         call. = FALSE)
  }
  if (class(input$signal) != "spec" | class(input$noise) != "spec") {
    stop("Input 'signal' and 'noise' must be of class 'spec'.", call. = FALSE)
  }
  if (length(input$signal$freq) != length(input$noise$freq)) {
    stop("'signal' and 'noise' spectra must have the same length.",
         call. = FALSE)
  }
  if (!all(input$signal$freq == input$noise$freq)) {
    stop("Frequency axes of 'signal' and 'noise' must match.")
  }
  if (!N >= 1) stop("'N' needs to be >= 1.", call. = FALSE)

  if (!is.null(freq.cut.lower)) {
    f1 <- which.min(abs(input$signal$freq - freq.cut.lower))
  }

  if (!is.null(freq.cut.upper)) {
    f2 <- which.min(abs(input$signal$freq - freq.cut.upper))
  } else {
    if (f2 == "max") {
      f2 <- length(input$signal$freq)
    }
  }

  freq <- input$signal$freq

  df <- diff(freq)
  df <- c(df[1], df)

  freq   <- freq[f1 : f2]
  df     <- df[f1 : f2]
  signal <- input$signal$spec[f1 : f2]
  noise  <- input$noise$spec[f1 : f2]

  snr <- list(
    freq = freq,
    spec = N * cumsum(signal * df) / cumsum(noise * df)
  )
  class(snr) <- "spec"

  return(snr)
}
