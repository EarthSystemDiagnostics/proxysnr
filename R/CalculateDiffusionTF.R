##' Diffusion transfer function
##'
##' This function implements an empirical Monte Carlo approach to estimate the
##' spectral transfer function for the effect of firn diffusion on the spatial
##' average of firn/ice-core stable isotope records.
##'
##' The approach is described in detail in Münch and Laepple (2018). In brief,
##' \code{nc} Gaussian white noise time series are created and diffused and the
##' average of these time series is calculated. The process is repeated
##' \code{ns} times. For each of the \code{ns} realisations, spectra of the
##' average diffused and undiffused records are calculated; subsequently, the
##' \code{ns} spectra are averaged, and the ratio of the average diffused to
##' the average undiffused spectrum yields the spectral transfer function.
##'
##' Diffusion is modelled as the convolution of the undiffused record with a
##' Gaussian with standard deviation given by the diffusion length
##' \code{sigma}. The spectral estimates are calculated using Thomson’s
##' multitaper method with three windows with linear detrending before
##' analysis.
##'
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
##' @param window length-2 vector giving a start and an end time (within `1 :
##' nt`) offering the possibility to only use a subset of the total length of
##' the simulated records for the transfer function analysis, while the default
##' of `NULL` means to use the records' entire lengths.
##' @param coherent if \code{TRUE}, \code{nc} identical white noise time series
##' are assumed to estimate the transfer function; else (the default) \code{nc}
##' independent noise series.
##' @param ... additional parameters which are passed to the spectral estimation
##' function \code{\link{SpecMTM}}.
##' @return a list of the components \code{signal}, \code{diffused} and
##' \code{ratio} which are objects of class \code{"spec"} providing the averages
##' over the \code{ns} simulations of:
##' \describe{
##' \item{\code{signal}:}{the undiffused noise spectrum}
##' \item{\code{diffused}:}{the diffused noise spectrum}
##' \item{\code{ratio}}{their ratio (diffused/undiffused), i.e. the transfer
##' function.}
##' }
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal- to centennial-scale isotope variations from Antarctic ice cores?
##' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
##' @export
CalculateDiffusionTF <- function(nt, nc, ns, sigma, res = 1, window = NULL,
                                 coherent = FALSE, ...) {

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
      if (window[1] < 1) stop("`window[1] must be >= 1.", call. = FALSE)
      if (window[2] > nt)
        stop("`window[2] is > total number of time points.", call. = FALSE)

      k <- window[1] : window[2]

    } else {

      k <- seq_len(nt)

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

            Xdiff[, j] <- DiffuseRecord(X[, j], sigma = sigma[, j], res = res)
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
    
    res <- list()
    res$signal     <- signal.spec.mean
    res$diffused   <- diff.spec.mean
    res$ratio$freq <- res$signal$freq
    res$ratio$spec <- res$diff$spec / res$signal$spec

    class(res$ratio) <- "spec"

    return(res)
    
}

