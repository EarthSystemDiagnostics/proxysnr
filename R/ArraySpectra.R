##' Spectral estimates of core array
##'
##' \code{ArraySpectra} calculates all relevant spectral estimates
##' from a given core array of \code{N} proxy records. The spectral estimates
##' can be smoothed in logarithmic space.
##'
##' The spectral estimates are calculated using Thomson’s multitaper method with
##' three windows with linear detrending before analysis
##' (see \code{\link{SpecMTM}}). Each spectral result is returned as an object
##' of class \code{"spec"} with the minimum elements \code{freq} and
##' \code{spec}.
##' @param cores a list of the proxy data from the core array. Each component is
##' expected to be a numeric vector of the proxy values of a common length.
##' @param res the sampling (e.g., temporal) resolution of the proxy data;
##' determines the frequency axis of the spectral estimates.
##' @param neff the effective number of records (e.g. to account for an expected
##' spatial correlation of the local noise). Per default, no spatial correlation
##' is assumed and \code{neff} is set to the number of proxy records (the length
##' of \code{cores}).
##' @param df.log width of the Gaussian kernel in logarithmic frequency units to
##' smooth the spectral estimates; \code{NULL} (the default) suppresses
##' smoothing.
##' @param ... additional parameters which are passed to the spectral estimation
##' function \code{\link{SpecMTM}}.
##' @return A list of the following components:
##' \describe{
##' \item{N:}{the number of (effective) proxy records of the core array.}
##' \item{single:}{a list of length \code{N} with the spectra of each individual
##' proxy record.}
##' \item{mean:}{the mean spectrum across all individual spectra.}
##' \item{stack:}{the spectrum of the average proxy record in the time domain
##' ("stacked record").}
##' }
##' @author Thomas Münch
##' @seealso \code{\link{SpecMTM}}
##' @export
ArraySpectra <- function(cores, res = 1, neff = length(cores),
                         df.log = NULL, ...) {

    # proxy data vectors must be of the same length
    if (sd(sapply(cores, function(lst) {length(lst)})) > 0) {
        stop("All data vectors in supplied input list must be of the same length.")
    }

    
    # estimate individual spectra
    single <- lapply(cores, function(lst) {
        SpecMTM(ts(lst, deltat = res), ...)
    })

    # estimate spectrum of stacked record
    ts.stack <- rowMeans(simplify2array(cores))
    stack <- SpecMTM(ts(ts.stack, deltat = res), ...)

    # calculate mean spectrum across individual record's spectra
    mean <- MeanSpectrum(single)


    # log-smooth spectra
    if (!is.null(df.log)) {

        single <- lapply(single, LogSmooth, df.log = df.log)
        mean   <- LogSmooth(mean, df.log = df.log)
        stack  <- LogSmooth(stack, df.log = df.log)
    }

    
    # return results as a list
    res <- list(
        N          = neff,
        single     = single,
        mean       = mean,
        stack      = stack)

    return(res)

}

