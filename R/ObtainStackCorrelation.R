#' Correlation with common signal
#'
#' This function calculates the theoretical correlation of a \code{"stacked"}
#' proxy record with the common signal depending on the number of records in
#' the stack and the time resolution of the records, given estimates of the
#' average proxy signal and noise spectra.
#'
#' The function is an implementation of Eqs. (6) and (7) in Münch and Laepple
#' (2018). The integral in (6) is approximated by the cumulative sum of the
#' integration arguments from \code{f.int1} to \code{f.int2}, where
#' \code{f.int1 = f1} and \code{f.int2} consecutively increases from \code{f1}
#' to \code{f2}.
#'
#' @param input a list of the spectral objects lists \code{signal} and
#'   \code{noise}, usually to be obtained from a call to
#'   \code{\link{SeparateSignalFromNoise}}.
#' @param N integer vector with the number of records in the assumed stack;
#'   correlations are then calculated for stacks with record numbers according
#'   to each element of \code{N}.
#' @param f1 index of the signal (and noise) frequency axis to specify the lower
#'   integration limit from which to integrate the spectra; per default the
#'   lowest frequency of the spectral estimates is omitted.
#' @param f2 as \code{f1} for the upper integration limit; defaults to use the
#'   maximum frequency of the given spectral estimates.
#' @param limits numeric vector with a frequency range of the integration: this
#'   is an alternative way of specifying the integration limits and overrides
#'   the setting by \code{f1} and \code{f2}. If not `NULL` it must be a length-2
#'   vector with the lower integration limit as first and the upper integration
#'   limit as second element.
#'
#' @return a list of two components:
#'   \describe{
#'   \item{freq:}{numeric vector of frequencies corresponding to the upper ends
#'     of the cumulative integrations;}
#'   \item{correlation:}{a \code{n * m} matrix where \code{n} corresponds to
#'     \code{length(N)} and \code{m} is given by \code{length(freq)} providing
#'     the correlation values as a function of the number of averaged records
#'     and the record resolution (= increasing upper frequency of the
#'     integration).}
#' }
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
ObtainStackCorrelation <- function(input, N = 1, f1 = 2, f2 = "max",
                                   limits = NULL) {

  snr <- GetIntegratedSNR(input, N = 1, f1 = f1, f2 = f2, limits = limits)

  getR <- function(n, x) {1 / sqrt(1 + 1 / (n * rev(x)))}

  correlation <- N %>%
    sapply(getR, x = snr$spec) %>%
    t()

  list(freq = snr$freq, correlation = correlation)

}

