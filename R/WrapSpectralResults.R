#' Wrap spectral results for several data sets
#'
#' This wrapper function is used to combine all main spectral results for the
#' signal, noise and signal-to-noise ratio as presented and discussed in Münch
#' and Laepple (2018), but it can also be used to combine the results for other
#' data sets.
#'
#' @param ... a comma separated list of named proxy data sets to analyse.
#' @param diffusion a list the same length as the number of datasets with each
#'   list element a spectral object of a transfer function (see also
#'   \code{\link{CalculateDiffusionTF}}) to correct the corresponding dataset
#'   for the effect of diffusion; i.e., internally the inverse of the transfer
#'   function values are used to correct for the diffusional smoothing (see
#'   Eq. 4 in Münch and Laepple, 2018). If \code{NULL}, no diffusion correction
#'   is applied at all; if instead you want to omit the diffusion correction
#'   only for some specific data set(s), set the corresponding list element(s)
#'   to \code{NA}.
#' @param time.uncertainty as \code{diffusion} a list of transfer functions to
#'   correct for the effect of time uncertainty
#'   (see also \code{\link{CalculateTimeUncertaintyTF}}).
#' @param df.log a vector of Gaussian kernel widths in log space to smooth the
#'   spectral estimates from each data set. If dimensions do not fit, its length
#'   is recycled to match the number of data sets.
#' @param crit.diffusion maximum diffusion correction value to obtain cutoff
#'   frequencies until which results are analysed to avoid large uncertainties
#'   at the high-frequency end of the spectra.
#'
#' @return A list of \code{N} lists, where \code{N} is the number of provided
#'   data sets and where each of these lists contains up to five elements:
#'   \describe{
#'   \item{\code{raw}:}{a list with four elements: three objects of class
#'     \code{"spec"} (the raw signal, noise and corresponding SNR spectra), and
#'     a two-element vector (\code{f.cutoff}) with the index and value of the
#'     cutoff frequency.}
#'   \item{\code{corr.diff.only}:}{as item \code{raw} but with the spectra after
#'     correction for the effect of diffusion.}
#'   \item{\code{corr.t.unc.only}:}{as item \code{raw} but with the spectra
#'     after correction for the effect of time uncertainty.}
#'   \item{\code{corr.full}:}{as item \code{raw} but with the spectra after
#'     correction for both the effects of diffusion and time uncertainty.}
#' }
#' The number of the returned list elements for each data set depends on
#' whether transfer functions for the corrections have been provided in
#' \code{diffusion} and \code{time.uncertainty} or not. Also, the element
#' \code{f.cutoff} is `NA` if diffusion has not been corrected for.
#'
#' @seealso \code{\link{CalculateDiffusionTF}},
#'   \code{\link{CalculateTimeUncertaintyTF}}
#' @author Thomas Münch
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @examples
#'
#' # Get main results of Münch and Laepple (2018)
#' 
#' results <- WrapSpectralResults(
#'                dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
#'                diffusion = diffusion.tf,
#'                time.uncertainty = time.uncertainty.tf,
#'                df.log = c(0.15, 0.15, 0.1))
#'
#' @export
#'
WrapSpectralResults <- function(..., diffusion = NULL, time.uncertainty = NULL,
                                df.log = 0.05, crit.diffusion = 2) {

  # Gather input data
  dat <- list(...)
  n <- length(dat)

  if (n == 0) stop("No data sets supplied.", call. = FALSE)

  if (is.null(diffusion)) diffusion <- as.list(rep(NA, n))
  if (is.null(time.uncertainty)) time.uncertainty <- as.list(rep(NA, n))

  # Check for correct dimensions
  if ((n != length(diffusion)) | (n != length(time.uncertainty))) {
    stop("Mismatch of dimensions of input data and correction function(s).")
  }
  
  if (length(df.log) == 1) df.log = rep(df.log, length.out = n)

  
  # Loop over data sets to obtain the relevant spectral quantities
  
  res <- list()
  for (i in 1 : n) {

    tmp <- list()

    # get diffusion/time-uncertainty correction functions
    d.flag <- !is.na(diffusion[i])
    if (d.flag) {
      d.crr <- 1 / diffusion[[i]]$spec
    }

    t.flag <- !is.na(time.uncertainty[i])
    if (t.flag) {            
      t.crr <- 1 / time.uncertainty[[i]]$spec
    }
    
    # critical cutoff frequency for diffusion correction
    if (d.flag) {
      idx <- which(d.crr >= crit.diffusion)[1]
      f.cutoff <- c(idx, diffusion[[i]]$freq[idx])
    }

    # mean and stack spectra
    spec <- ObtainArraySpectra(dat[[i]], df.log = df.log[i])

    # raw signal and noise spectra
    tmp$raw <- SeparateSignalFromNoise(spec)
    tmp$raw$f.cutoff <- NA_real_

    # corrected signal and noise spectra
    if (d.flag) {
      tmp$corr.diff.only <-
        SeparateSignalFromNoise(spec, diffusion = d.crr)
      tmp$corr.diff.only$f.cutoff <- f.cutoff
    }
    if (t.flag) {
      tmp$corr.t.unc.only <-
        SeparateSignalFromNoise(spec, time.uncertainty = t.crr)
      tmp$corr.t.unc.only$f.cutoff <- NA_real_
    }
    if (d.flag & t.flag) {
      tmp$corr.full <-
        SeparateSignalFromNoise(spec, time.uncertainty = t.crr,
                                diffusion = d.crr)
      tmp$corr.full$f.cutoff <- f.cutoff
    }

    res[[i]] <- tmp

  }

  
  # Output
  
  names(res) <- names(dat)

  return(res)

}

