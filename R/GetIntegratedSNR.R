#' Frequency-dependent signal-to-noise ratio from integrated spectra
#'
#' This function calculates the signal-to-noise ratio (SNR) as a function of
#' frequency, interpreted as the temporal resolution of a proxy record. This
#' variant of the SNR is obtained from a signal and a noise spectrum which are
#' each integrated across frequencies before taking their ratio.
#'
#' The function is an implementation of Eq. (6) in Münch and Laepple
#' (2018). The integral in (6) is approximated by the cumulative sum of the
#' integration arguments from \code{f.int1} to \code{f.int2}, where
#' \code{f.int1 = f1} and \code{f.int2} consecutively increases from \code{f1}
#' to \code{f2}.
#'
#' @param input a list of the spectral objects lists \code{signal} and
#'   \code{noise}, usually to be obtained from a call to
#'   \code{\link{SeparateSignalFromNoise}}.
#' @param N integer; number of proxy records averaged. The default returns the
#'   SNR of a single proxy record. For a different number, the SNR is calculated
#'   for a "stack" averaged across \code{N} individual proxy records with the same
#'   signal and equivalent noise characteristics, assuming independent noise
#'   between the records.
#' @inheritParams ObtainStackCorrelation
#'
#' @return a spectral object list of the SNR.
#'
#' @author Thomas Münch
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @export
#'
GetIntegratedSNR <- function(input, N = 1, f1 = 2, f2 = "max", limits = NULL) {

  # error checking

  if (!is.list(input)) stop("`input` must be a list.", call. = FALSE)

  if (!all(utils::hasName(input, c("signal", "noise"))))
    stop("`input` must have elements `signal` and `noise`.", call. = FALSE)

  is.spectrum(input$signal)
  is.spectrum(input$noise)

  if (length(input$signal$freq) != length(input$noise$freq)) {
    stop("`signal` and `noise` must have the same number of spectral estimates.",
         call. = FALSE)
  }

  if (!all(input$signal$freq == input$noise$freq)) {
    stop("Frequency axes of `signal` and `noise` do not match.")
  }

  if (length(N) != 1) stop("`N` must have length 1.", call. = FALSE)
  if (N < 1) stop("`N` must be >= 1.", call. = FALSE)

  if (!is.null(limits)) {

    if (length(limits) != 2) stop("`limits` must have length 2.", call. = FALSE)
    if (diff(limits) <= 0)
      stop("`limits[2]` must be > `limits[1]`.", call. = FALSE)

    f1 <- which.min(abs(input$signal$freq - limits[1]))
    f2 <- which.min(abs(input$signal$freq - limits[2]))

  } else {

    if (length(f1) != 1) stop("`f1` must have length 1.", call. = FALSE)
    if (length(f2) != 1) stop("`f2` must have length 1.", call. = FALSE)

    if (f1 < 1) stop("`f1` must be >= 1.", call. = FALSE)

    if (is.character(f2)) {
      if (f2 != "max")
        stop("`f2` must be set to \"max\" or to an index.", call. = FALSE)
      f2 <- length(input$signal$freq)
    } else {
      if (f2 > length(input$signal$freq))
        stop("`f2` larger than length of data.", call. = FALSE)
    }

    if (f2 <= f1) stop("`f2` must be > `f1`.", call. = FALSE)

  }

  # do the integration

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
