##' Final DML and WAIS spectra
##'
##' Internal function to produce the final DML and WAIS spectra (signal, noise,
##' SNR) used for Figs. (3) and (4) in Münch and Laepple (2018). The DML spectra
##' are a combination of the results from the DML1 and DML2 data sets.
##' @param spec output from \code{\link{WrapSpectralResults}} for DML and WAIS
##' data.
##' @param dml.knit.f frequency at which to combine the spectra from the DML1
##' and DML2 data sets (defaults to 1/10 yr^(-1)); DML2 spectra are used for
##' lower frequencies than \code{dml.knit.f}, DML1 for the higher frequencies.
##' @param df.log width of Gaussian kernel in logarithmic frequency units for
##' additional smoothing for visual purposes; \code{NULL} suppresses smoothing.
##' @return A list with the components \code{dml} and \code{wais}, where each is
##' a list containing three elements of class \code{"spec"} and one numeric
##' vector:
##' \describe{
##' \item{\code{signal}:}{the signal spectrum.}
##' \item{\code{noise}:}{the noise spectrum.}
##' \item{\code{snr}:}{the frequency-dependent signal-to-noise ratio.}
##' \item{\code{f.cutoff}:}{a two-element vector with the index and value of the
##' cutoff frequency from constraining the diffusion correction.}
##' }
##' @author Thomas Münch
##' @references Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal- to centennial-scale isotope variations from Antarctic ice cores?
##' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
PublicationSNR <- function(spec, dml.knit.f = 0.1, df.log = 0.125) {

    
    # Combine DML1 and DML2 spectra
    
    idx.knit1 <- which(spec$dml1$raw$signal$freq > dml.knit.f)
    idx.knit2 <- which(spec$dml2$raw$signal$freq < dml.knit.f)
    
    dml.signal <- list()
    dml.signal$freq <- c(spec$dml2$raw$signal$freq[idx.knit2],
                         spec$dml1$raw$signal$freq[idx.knit1])
    dml.signal$spec <- c(spec$dml2$corr.full$signal$spec[idx.knit2],
                         spec$dml1$corr.full$signal$spec[idx.knit1])

    dml.noise <- dml.signal
    dml.noise$spec <- c(spec$dml2$corr.full$noise$spec[idx.knit2],
                        spec$dml1$corr.full$noise$spec[idx.knit1])

    dml.f.cutoff <- spec$dml1$f.cutoff
    dml.f.cutoff[1] <- which(dml.signal$freq == spec$dml1$f.cutoff[2])


    # Smooth DML & WAIS signal and noise and calculate SNR

    if (!is.null(df.log)) {
        
        dml.signal <- LogSmooth(dml.signal, df.log = df.log)
        dml.noise  <- LogSmooth(dml.noise, df.log = df.log)
    
        wais.signal <- LogSmooth(spec$wais$corr.full$signal, df.log = df.log)
        wais.noise  <- LogSmooth(spec$wais$corr.full$noise, df.log = df.log)
    }

    dml.snr <- list()
    dml.snr$freq <- dml.signal$freq
    dml.snr$spec <- dml.signal$spec / dml.noise$spec

    wais.snr <- list()
    wais.snr$freq <- wais.signal$freq
    wais.snr$spec <- wais.signal$spec / wais.noise$spec

    
    # Organize output

    res <- list()

    class(dml.signal) <- "spec"
    class(dml.noise)  <- "spec"
    class(dml.snr)    <- "spec"

    class(wais.signal) <- "spec"
    class(wais.noise)  <- "spec"
    class(wais.snr)    <- "spec"

    res$dml$signal   <- dml.signal
    res$dml$noise    <- dml.noise
    res$dml$snr      <- dml.snr
    res$dml$f.cutoff <- dml.f.cutoff
    
    res$wais$signal   <- wais.signal
    res$wais$noise    <- wais.noise
    res$wais$snr      <- wais.snr
    res$wais$f.cutoff <- spec$wais$f.cutoff

    return(res)

}

