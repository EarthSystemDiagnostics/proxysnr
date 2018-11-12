##' Diffusion transfer function
##'
##' This function implements an empirical Monte Carlo approach to estimate the
##' spectral transfer function for the effect of firn diffusion on the spatial
##' average of isotope records.
##'
##' The approach is described in detail in Münch and Laepple (2018). In brief,
##' \code{nc} white noise time series are created and diffused and the average
##' of these time series is calculated. The process is repeated \code{ns}
##' times. For each of the \code{ns} realisations, spectra of the average
##' diffused and undiffused records are calculated yielding the spectral
##' transfer function.
##'
##' Diffusion is modelled as the convolution of the undiffused record with a
##' Gaussian with standard deviation given by the diffusion length. The
##' spectral estimates are calculated using the function \code{SpecMTM} from the
##' package \code{PaleoSpec}, which applies the multi-taper method. If the
##' package is not available, spectra are calculated using base \code{R}
##' functions.
##' @param nt the length of the modelled isotope records (i.e. the number of
##' data points in each record)
##' @param nc the number of cores in the modelled core array
##' @param ns the number of Monte Carlo simulations for estimating the transfer
##' function
##' @param sigma numeric vector of length \code{nt} or numeric array of
##' dimension \code{nt * nc} providing diffusion length values. The \code{nt}
##' diffusion length values are assumed to correspond to the respective
##' \code{nt} isotope values. If only a numeric vector is provided, it is
##' assumed that these diffusion lengths are valid for all \code{nc} cores. If
##' an array is provided, each column provides the diffusion lengths for the
##' respective core. Note that the units of the diffusion length must match the
##' units of \code{res}.
##' @param res the sampling (e.g., temporal) resolution of the isotope data;
##' determines the frequency axis of the transfer function.
##' @param coherent if \code{TRUE}, \code{nc} identical white noise time series
##' are assumed to estimate the transfer function; else (the default) \code{nc}
##' independent noise series.
##' @return a list of the components \code{signal}, \code{diffused} and
##' \code{ratio} which are \code{"spectral object"} lists providing averages
##' over the \code{ns} simulations of:
##' \describe{
##' \item{\code{signal}:}{the undiffused noise spectrum}
##' \item{\code{diffused}:}{the diffused noise spectrum}
##' \item{\code{ratio}}{their ratio (diffused/undiffused), i.e. the transfer
##' function.}
##' }
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal to centennial scale isotope variations from Antarctic ice cores?
##' Clim. Past Discuss., https://doi.org/10.5194/cp-2018-112, in review, 2018.
##' @export
DiffusionTF <- function(nt, nc, ns, sigma, res = 1, coherent = FALSE) {

    # convert sigma vector into array if necessary and check for dimensions
    if (is.null(dim(sigma))) {
        
        if (length(sigma != nt)) {
            stop("Invalid length of supplied vector of diffusion lengths.")
        }
        
        sigma <- array(sigma, dim = c(nt, nc))

    } else {

        if (prod(dim(sigma)) != (nt * nc)) {
            stop("Invalid dimensions of supplied array of diffusion lengths.")
        }

    }

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


    # create 'nc' white noise time series, diffuse them, and calculate their
    # average; repeat 'ns' times
    
    stack.signal <- array(dim = c(nt, ns))
    stack.diff <- array(dim = c(nt, ns))

    for (i in 1 : ns) {

        if (coherent) {
            X <- array(rnorm(nt), dim = c(nt, nc))
        } else {
            X <- array(rnorm(nt * nc), dim = c(nt, nc))
        }
        
        Xdiff <- array(dim = c(nt, nc))

        for (j in 1 : nc) {

            Xdiff[, j] <- DiffuseRecord(X[, j], sigma = sigma[, j], res = res)
        }

        stack.signal[, i] <- rowMeans(X)
        stack.diff[, i] <- rowMeans(Xdiff)

    }


    # calculate spectra of the 'ns' average diffused and undiffused noise series

    signal.lst <- lapply(seq_len(ncol(stack.signal)), function(i) {
        stack.signal[, i]})
    signal.spec <- lapply(signal.lst, function(x) {
        fun(ts(x, deltat = res), ...)})

    diff.lst <- lapply(seq_len(ncol(stack.diff)), function(i) {
        stack.diff[, i]})
    diff.spec <- lapply(diff.lst, function(x) {
        fun(ts(x, deltat = res), ...)})


    # calculate the average over all 'ns' simulations

    if (paleospec) {
        signal.spec.mean <- PaleoSpec::MeanSpectrum(signal.spec, iRem = 0)$spec
        diff.spec.mean  <- PaleoSpec::MeanSpectrum(diff.spec, iRem = 0)$spec
    } else {
        signal.spec.mean <- list()
        signal.spec.mean$freq <- signal.spec[[1]]$freq
        signal.spec.mean$spec <- rowMeans(sapply(signal.spec,
                                                 function(lst) {lst$spec}))
        diff.spec.mean <- list()
        diff.spec.mean$freq <- diff.spec[[1]]$freq
        diff.spec.mean$spec <- rowMeans(sapply(diff.spec,
                                               function(lst) {lst$spec}))
        
    }


    # return average undiffused and diffused spectra and the corresponding
    # ratio (transfer function) as a list
    
    res <- list()
    res$signal     <- signal.spec.mean
    res$diffused   <- diff.spec.mean
    res$ratio$freq <- res$signal$freq
    res$ratio$spec <- res$diff$spec / res$signal$spec

    return(res)
    
}

