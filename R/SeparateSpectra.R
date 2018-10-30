##' Calculate signal and noise spectra
##' 
##' \code{SeparateSpectra} calculates the raw signal and noise spectra, and the
##' corresponding signal-to-noise ratio, as estimated from a core array of
##' \code{N} proxy records, and allows one to correct these, where applicable,
##' for the effects of time uncertainty and water vapour diffusion (relevant for
##' firn and ice cores).
##'
##' This function is an implementation of Eq. (4) in Münch and Laepple
##' (2018). While the diffusion correction is relevant only for diffusing
##' proxies such as stable isotopes from firn and ice cores, this function can
##' be applied to a large set of proxy data since only one of the two
##' correction functions, or none, need to be supplied; thus, e.g., it can also
##' be used for non-diffusing proxy data where only time uncertainty is
##' relevant, or for estimating raw signal and noise spectra by supplying no
##' correction functions at all.
##' @param spectra a list of the raw spectral estimates from a proxy core
##' array. Expected is the output from \code{\link{ArraySpectra}}, but
##' sufficient is a named list of two components giving the \code{mean} and
##' \code{stack} spectra.
##' @param neff the effective number of records (e.g. to account for an expected
##' spatial correlation of the local noise). Per default set to element \code{N}
##' in \code{spectra}, otherwise supply it explicitly here.
##' @param diffusion numeric vector of diffusion correction values (inverse
##' transfer function); must be of the same length as the spectral estimates in
##' \code{spectra}. By omitting this parameter no correction will be applied.
##' @param time.uncertainty numeric vector of time uncertainty correction
##' values (inverse transfer function); must be of the same length as the
##' spectral estimates in \code{spectra}. By omitting this parameter no
##' correction will be applied.
##' @return A list of three components:
##' \describe{
##' \item{\code{signal}:}{the raw or corrected signal spectrum.}
##' \item{\code{noise}:}{the raw or corrected noise spectrum.}
##' \item{\code{snr}:}{the signal-to-noise ratio as calculated from the previous
##' components.}
##' }
##' @author Thomas Münch
##' @seealso \code{\link{ArraySpectra}}
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal to centennial scale isotope variations from Antarctic ice cores?
##' Clim. Past Discuss., https://doi.org/10.5194/cp-2018-112, in review, 2018.
##' @export
SeparateSpectra <- function(spectra, neff = spectra$N,
                            diffusion, time.uncertainty) {

    # error checking

    if (is.null(neff)) {
        stop("Supply (effective) number of records.")
    }

    if (missing(diffusion)) {
        diffusion <- NULL
    }

    if (missing(time.uncertainty)) {
        time.uncertainty <- NULL
    }

    if (any(is.na(match(c("mean", "stack"), names(spectra))))) {
        stop("Input list of spectra does not seem to have the right format.")
    }

    if (is.null(diffusion)) {
        diffusion <- 1
    } else {
        if (length(diffusion) != length(spectra$mean$freq)) {
            stop("Length of diffusion correction does not match length of
        spectral estimates.")
        }
    }

    if (is.null(time.uncertainty)) {
        time.uncertainty <- 1
    } else {
        if (length(diffusion) != length(spectra$mean$freq)) {
            stop("Length of time uncertainty correction does not match length of
        spectral estimates.")
        }
    }


    # calculate the signal and noise spectra
    
    d.corr <- diffusion
    t.corr <- time.uncertainty

    N <- neff
    n <- N / (N - t.corr)

    corr.fac <- n * d.corr

    mean  <- spectra$mean
    stack <- spectra$stack
    
    signal <- list()
    noise  <- list()
    snr    <- list()

    signal$freq <- mean$freq
    noise$freq  <- mean$freq
    snr$freq    <- mean$freq

    signal$spec <- t.corr * corr.fac * (stack$spec - mean$spec / N)
    noise$spec  <- corr.fac * (mean$spec - t.corr * stack$spec)

    snr$spec <- signal$spec / noise$spec

    
    # return results as a list
    res <- list(
        signal  = signal,
        noise   = noise,
        snr     = snr)

    return(res)

}

