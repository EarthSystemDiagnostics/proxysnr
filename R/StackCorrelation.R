##' Correlation with common signal
##'
##' This function calculates the theoretical correlation of a \code{"stacked"}
##' proxy record with the common signal depending on the number of records in
##' the stack and the time resolution of the records, given estimates of the
##' average proxy signal and noise spectra.
##'
##' The function is an implementation of Eqs. (6) and (7) in Münch and Laepple
##' (2018). The integral in (6) is approximated by the cumulative sum of the
##' integration arguments from \code{f.int1} to \code{f.int2}, where
##' \code{f.int1 = f1} and \code{f.int2} consecutively increases from \code{f1}
##' to \code{f2}.
##' @param input a list of the spectral objects lists \code{signal} and
##' \code{noise}, usually to be obtained from a call to
##' \code{\link{SeparateSpectra}}
##' @param N maximum number of records in the assumed stack; correlations are
##' then calculated for stacks with record numbers from 1 to \code{N}
##' @param f1 index of the the minimum frequency from which to integrate the
##' signal and noise spectra for calculating the correlation; per default the
##' lowest frequency of the spectral estimates is omitted
##' @param f2 index of the maximum frequency until which to integrate the signal
##' and noise spectra for calculating the correlation; defaults to use the
##' maximum frequency of the given spectral estimates
##' @param freq.cut.lower lower frequency (not index!) at which to cut the
##' spectra: this provides a direct way for specifying a minimum frequency for
##' the integration different from the minimum frequency of the spectral
##' estimates. Setting \code{freq.cut.lower} overrides the frequency
##' corresponding to the index set in \code{f1}.
##' @param freq.cut.upper upper frequency (not index!) at which to cut the
##' spectra: this provides a direct way for specifying a maximum frequency for
##' the integration different from the maximum frequency of the spectral
##' estimates. Setting \code{freq.cut.upper} overrides the frequency
##' corresponding to the index set in \code{f2}.
##' @return a list of three components:
##' \describe{
##' \item{signal:}{a spectral object list of the integrated signal}
##' \item{noise:}{a spectral object list of the integrated noise}
##' \item{correlation:}{a \code{n * m} matrix where \code{n} corresponds to
##' \code{N} and \code{m} is given by \code{length(f1 : f2)} providing the
##' correlation values as a function of the number of averaged records and
##' the record resolution (= increasing upper frequency of the integration)}
##' }
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal- to centennial-scale isotope variations from Antarctic ice cores?
##' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
##' @export
StackCorrelation <- function(input, N, f1 = 2, f2 = "max",
                             freq.cut.lower = NULL, freq.cut.upper = NULL) {

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

    df <- diff(input$signal$freq)
    df <- c(df[1], df)

    input$signal$freq <- input$signal$freq[f1 : f2]
    input$signal$spec <- input$signal$spec[f1 : f2]
    input$noise$freq  <- input$noise$freq[f1 : f2]
    input$noise$spec  <- input$noise$spec[f1 : f2]

    df <- df[f1 : f2]

    signal.int <- noise.int <- list()
    signal.int$freq <- noise.int$freq <- input$signal$freq
    signal.int$spec <- 2 * cumsum(input$signal$spec * df)
    noise.int$spec  <- 2 * cumsum(input$noise$spec * df)

    snr.int <- signal.int$spec / noise.int$spec

    correlation <- matrix(nrow = N, ncol = length(snr.int))
    for (i in 1 : N) {
        
        correlation[i, ] <- 1 / sqrt(1 + 1 / (i * rev(snr.int)))
    }

    res <- list()
    res$signal <- signal.int
    res$noise  <- noise.int
    res$correlation <- correlation
    
    return(res)

}
        
