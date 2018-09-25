##' Spectral estimates of core array
##'
##' \code{ArraySpectra} calculates all relevant spectral estimates
##' from a given core array of \code{N} proxy records. All spectral estimates
##' can be smoothed in logarithmic space and confidence intervals can be
##' calculated.
##'
##' The spectral estimates are calculated using the function
##' \code{SpecMTM} from the package \code{PaleoSpec}, which applies the
##' multi-taper method. If the package is not available, spectra are
##' calculated using base \code{R} functions, which may result in estimates with
##' higher spectral uncertainty. In both cases, each spectral result is returned
##' as a \code{"spectral object"} list with the minimum elements \code{freq} and
##' \code{spec}. Log-smoothing and (re-)calculation of confidence intervals is
##' only available with \code{PaleoSpec}.
##'
##' \code{PaleoSpec} is available upon request.
##' @param cores a list of the proxy data from the core array. Each component is
##' expected to be a numeric vector of the proxy values of a common length.
##' @param res the sampling (e.g., temporal) resolution of the proxy data;
##' determines the frequency axis of the spectral estimates.
##' @param neff the effective number of records to account for an expected
##' spatial correlation of the local noise. Per default set to the number of
##' proxy records (the length of \code{cores}).
##' @param df.log width of the Gaussian kernel in logarithmic
##' space to smooth the spectral estimates; \code{NULL} (the default) suppresses
##' smoothing. Log-smoothing automatically results in calculation of
##' confidence intervals using a p-value of \code{0.05}.
##' @param pval p-value for calculating confidence intervals of the spectral
##' estimates if no log-smoothing is applied to the spectra, or for
##' re-calculating the confidence intervals with a different p-value than the
##' default \code{0.05}.
##' @param ... additional parameters supplied to the spectral estimation
##' function.
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
##' @seealso \code{[PaleoSpec]{SpecMTM}}
##' @export
ArraySpectra <- function(cores, res = 1, neff = length(cores),
                         df.log = NULL, pval = 0.05, ...) {

    # check if package PaleoSpec is available
    if (!requireNamespace("PaleoSpec", quietly = TRUE)) {
        paleospec <- FALSE
        fun <- spectrum
        warning(paste("Package \"PaleoSpec\" not found.",
                      "Using base \"R\" methods for spectral estimation."),
                call. = FALSE)
    } else {
        paleospec <- TRUE
        fun <- PaleoSpec::SpecMTM
    }

    # proxy data vectors must be of the same length
    if (sd(sapply(cores, function(lst) {length(lst)})) > 0) {
        stop("All data vectors in supplied input list must be of the same length.")
    }

    # estimate individual spectra
    single <- lapply(cores, function(lst) {
        fun(ts(lst, deltat = res), ...)
    })

    # estimate spectrum of stacked record
    ts.stack <- rowMeans(simplify2array(cores))
    stack <- fun(ts(ts.stack, deltat = res), ...)

    # calculate mean spectrum across individual record's spectra
    if (paleospec) {
        mean <- PaleoSpec::MeanSpectrum(single, iRemoveLowest = 0)$spec
    } else {
        mean <- list()
        mean$freq <- stack$freq
        mean$spec <- rowMeans(sapply(single, function(lst) {lst$spec}))
    }


    # log-smooth spectra
    if (!is.null(df.log)) {

        if (paleospec) {

            single <- lapply(single, PaleoSpec::LogSmooth, df.log = df.log)
            mean   <- PaleoSpec::LogSmooth(mean, df.log = df.log)
            stack  <- PaleoSpec::LogSmooth(stack, df.log = df.log)
            
        } else {

            warning("Log-smoothing not possible without package \"PaleoSpec\".")

        }
    }

    
    # re-set confidence intervals
    if (is.null(df.log) | (!is.null(df.log) & pval != 0.05)) {

        if (paleospec) {
            
            single <- lapply(single, function(lst) {
                lst <- PaleoSpec::AddConfInterval(lst, pval = pval)})
            mean  <- PaleoSpec::AddConfInterval(mean, pval = pval)
            stack <- PaleoSpec::AddConfInterval(stack, pval = pval)

        } else {

            warning(paste("(Re-)calculation of confidence intervals",
                          "not possible without package \"PaleoSpec\"."))
        }
    }

    
    # return results as a list
    res <- list(
        N          = N,
        single     = single,
        mean       = mean,
        stack      = stack)

    return(res)

}

