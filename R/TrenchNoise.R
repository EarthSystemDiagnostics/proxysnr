##' Trench noise spectrum
##'
##' Internal function to obtain estimates of the noise spectrum at the local
##' scale at EDML using the available isotope data from the shallow trench
##' profiles from Münch et al. (2017). To convert the trench depth axis into a
##' time axis, a constant accumulation rate with a given uncertainty is
##' assumed.
##' @param trench.data list of the shallow trench isotope profiles from Münch et
##' al. (2017)
##' @param acc.rate assumed mean accumulation rate at the trench site [cm/yr]
##' @param sigma.acc.rate assumed uncertainty of the accumulation rate [cm/yr]
##' @param res depth resolution of the trench data; defaults to 3 cm.
##' @param neff effective number of trench records to account for the spatial
##' autocorrelation of the stratigraphic noise; defaults to 19 for the total
##' available number of 22 profiles (Münch et al., 2017).
##' @param df.log Gaussian kernel width in log space to smooth the spectral
##' estimate; defaults to 0.1.
##' @return A list of the three components:
##' \describe{
##' \item{\code{noise.lower}:}{spectral object list with the estimated trench
##' noise spectrum for the assumed lower bound of the accumulation rate.}
##' \item{\code{noise.mean}:}{as previous but for the assumed mean accumulation
##' rate.}
##' \item{\code{noise.upper}:}{as previous but for the assumed upper bound of
##' the accumulation rate.}
##' }
##' @author Thomas Münch
##' @references
##' Münch, T. and Laepple, T.: What climate signal is contained in
##' decadal to centennial scale isotope variations from Antarctic ice cores?
##' Clim. Past Discuss., https://doi.org/10.5194/cp-2018-112, in review, 2018.
TrenchNoise <- function(trench.data, acc.rate = 25, sigma.acc.rate = 5,
                        res = 3, neff = 19, df.log = 0.1) {

    
    # Calculate raw (w/o diffusion correction) noise spectra depending on
    # accumulation rate
    
    res <- list()

    # calculation for lower bound of accumulation rate
    res$noise.lower <- SeparateSpectra(
        ArraySpectra(
            trench.data,
            res = res / (acc.rate - sigma.acc.rate),
            neff = neff,
            df.log = df.log
        )
    )
    # calculation for mean of accumulation rate
    res$noise.mean <- SeparateSpectra(
        ArraySpectra(
            trench.data,
            res = res / (acc.rate),
            neff = neff,
            df.log = df.log
        )
    )
    # calculation for upper bound of accumulation rate
    res$noise.lower <- SeparateSpectra(
        ArraySpectra(
            trench.data,
            res = res / (acc.rate + sigma.acc.rate),
            neff = neff,
            df.log = df.log
        )
    )

    return(res)

}

