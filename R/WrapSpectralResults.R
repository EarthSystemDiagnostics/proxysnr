##' Main results for Münch and Laepple (2018)
##'
##' This internal wrapper function is used to combine all main spectral results
##' presented and discussed in Münch and Laepple (2018) but it can also be used
##' to combine the results for other data sets.
##' @param ... a comma separated list of named proxy data sets to analyse.
##' @param diffusion a list of (inverse) transfer functions to correct for the
##' effect of diffusion. The length of the list has to match the number of
##' provided data sets, thus, one transfer function per data set is assumed. 
##' @param time.uncertainty similar to \code{diffusion} a list of (inverse)
##' transfer functions to correct for the effect of time uncertainty.
##' @param df.log a vector of Gaussian kernel widths in log space to smooth the
##' spectral estimates from each data set. If dimensions do not fit, its length
##' is recycled to match the number of data sets.
##' @param crit.diffusion maximum diffusion correction value to obtain cutoff
##' frequencies until which results are analysed to avoid large uncertainties at
##' the high-frequency end of the spectra.
##' @param inverse.tf logical; if \code{TRUE}, it is assumed that
##' \code{diffusion} and \code{time.uncertainty} provide the inverse transfer
##' functions which can be readily used to correct the spectra. If \code{FALSE}
##' (the default), the inverse of the provided transfer functions is calculated
##' within the function and used for the corrections. See Eqs. (4) in Münch and
##' Laepple (2018) for the definitions.
##' @return A list of \code{N} lists, where \code{N} is the number of provided
##' data sets and where each of these lists contains four elements:
##' \describe{
##' \item{\code{f.cutoff}:}{a two-element vector with the index and value of the
##' cutoff frequency.}
##' \item{\code{raw}:}{the raw signal, noise and corresponding SNR spectra.}
##' \item{\code{corr.t.unc.only}:}{the signal, noise and corresponding SNR
##' spectra after correction for the effect of time uncertainty.}
##' \item{\code{corr.full}:}{the signal, noise and corresponding SRN spectra
##' after correction for both the effects of diffusion and time uncertainty.}
##' }
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal to centennial scale isotope variations from Antarctic ice cores?
##' Clim. Past Discuss., https://doi.org/10.5194/cp-2018-112, in review, 2018.
##' @examples
##' # Get main results of Münch and Laepple (2018)
##' 
##' results <- proxysnr:::WrapSpectralResults(
##'                dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
##'                diffusion = diffusion.tf,
##'                time.uncertainty = time.uncertainty.tf,
##'                df.log = c(0.15, 0.15, 0.1))
WrapSpectralResults <- function(..., diffusion, time.uncertainty,
                                df.log = 0.05, crit.diffusion = 2,
                                inverse.tf = FALSE) {

    # Gather input data
    dat <- list(...)

    # Check for correct dimensions
    n <- length(dat)
    if ((n != length(diffusion)) | (n != length(time.uncertainty))) {
        stop("Mismatch of dimensions of input data and corretion function(s).")
    }
    
    if (length(df.log) == 1) df.log = rep(df.log, length.out = n)

    
    # Loop over data sets to obtain the relevant spectral quantities
    
    res <- list()
    for (i in 1 : n) {

        tmp <- list()

        # get diffusion/time-uncertainty correction functions
        d.crr <- diffusion[[i]]$spec
        t.crr <- time.uncertainty[[i]]$spec
        if (!inverse.tf) {
            d.crr <- 1 / d.crr
            t.crr <- 1 / t.crr
        }

        # critical cutoff frequency for diffusion correction
        idx <- which(d.crr >= crit.diffusion)[1]
        tmp$f.cutoff <- c(idx, diffusion[[i]]$freq[idx])

        # mean and stack spectra
        spec <- ArraySpectra(dat[[i]], df.log = df.log[i])

        # raw signal and noise spectra
        tmp$raw <- SeparateSpectra(spec)

        # corrected signal and noise spectra
        tmp$corr.t.unc.only <-
            SeparateSpectra(spec, time = t.crr)
        tmp$corr.full <-
            SeparateSpectra(spec, time = t.crr, diffusion = d.crr)

        res[[i]] <- tmp

    }

    
    # Output
    
    names(res) <- names(dat)

    return(res)

}

