#' DML oxygen isotope data
#'
#' A dataset containing the oxygen isotope data of 15 firn cores from Dronning
#' Maud Land (DML), Antarctica, as used in Münch and Laepple (2018).
#'
#' @format A named list with two components:
#' \describe{
#'   \item{dml1:}{a named list with 15 components where each component is a
#'     vector of length 194 with the oxygen isotope data in per mil for all 15
#'     cores for the years from 1994 CE to 1801 CE.}
#'   \item{dml2:}{a named list with 3 components where each component is a vector
#'     of length 995 with the oxygen isotope data in per mil for the cores
#'     \code{"B31"}, \code{"B32"} and \code{"B33"} for the years from 1994 CE to
#'     1000 CE.}
#' }
#' @details
#' Data has been processed from the original source data as documented in the
#' package source (\code{data-raw/load.DML.R}). Cores included are FB9804,
#' FB9805, FB9807–FB9811, FB9813–FB9817, and B31–B33.
#' @source
#' Data citation:\cr
#' Graf, W., et al.: Stable-isotope records from Dronning Maud Land, Antarctica,
#' \url{https://doi.org/10.1594/PANGAEA.728240}, 2002.
#'
#' Literature citation:\cr
#' Graf, W., et al.: Stable-isotope records from Dronning Maud Land, Antarctica,
#' Ann. Glaciol., 35, 195–201, 2002.\cr
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#' 2018.
"dml"

#' WAIS oxygen isotope data
#'
#' A dataset containing the oxygen isotope data of 5 firn cores from the West
#' Antarctic Ice Sheet (WAIS), Antarctica, as used in Münch and Laepple (2018).
#'
#' @format A data frame with 201 rows and 5 variables where each variable is the
#' time series of oxygen isotope data from one of the firn cores in per mil for
#' the years from 2000 CE to 1800 CE.
#' @details
#' Data has been processed from the original source data as documented in the
#' package source (\code{data-raw/load.WAIS.R}). Cores included are WDC2005A,
#' ITASE-1999-1, ITASE-2000-1, ITASE-2000-4 and ITASE-2000-5.
#' @source
#' Data citation:\cr
#' Steig, E. J.: West Antarctica Ice Core and Climate Model Data. WDC2005A,
#' ITASE-1999-1, ITASE-2000-1, ITASE-2000-4, ITASE-2000-5. Boulder, Colorado
#' USA: National Snow and Ice Data Center,
#' \url{https://doi.org/10.7265/N5QJ7F8B}, 2013.
#'
#' Literature citation:\cr
#' Steig, E. J., et al.: Recent climate and ice-sheet changes in West Antarctica
#' compared with the past 2,000 years, Nat. Geosci., 6, 372–375,
#' \url{https://doi.org/10.1038/ngeo1778}, 2013.\cr
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#' 2018.
"wais"

#' Diffusion transfer functions
#'
#' A dataset containing the calculated diffusion transfer functions for the
#' DML1, DML2 and WAIS records from Münch and Laepple (2018) based on numerical
#' simulations. The transfer function describes the average effect firn
#' diffusion has on the mean of the spectra from a number of `N` spatially
#' distributed oxygen isotope records.
#'
#' @format A named list with the 3 elements \code{dml1}, \code{dml2}, and
#'   \code{wais}, where each element is a spectral object:
#'   \describe{
#'   \item{\code{freq}:}{the frequency axis in units of \code{1 / yr},}
#'   \item{\code{spec}:}{the transfer function value at each frequency.}
#' }
#'
#' @source
#' The transfer functions were obtained using \code{?CalculateDiffusionTF} for
#' the site-specific diffusion lengths provided by
#' \code{proxysnr:::diffusion.length}; see also the respective package
#' vignette \code{vignette("calculate-transfer-functions")}. Isotope records
#' were simulated at semi-annual resolution and the transfer functions
#' interpolated in frequency space to annual resolution. Note here that prior to
#' `proxysnr` v1.0.0, the transer functions supplied to
#' `SeparateSignalFromNoise()` needed to match in frequency axis the spectra of 
#' the analysed data, hence the interpolation step. From `proxysnr` v1.0.0
#' onwards this is no longer required, as the frequency axis overlap is checked
#' within `SeparateSignalFromNoise()`, with interpolation taking place there, if
#' needed. The example firn diffusion transfer functions here were, however,
#' left as originally produced to maintain consistency with the results
#' published in Münch and Laepple (2018).
#'
#' @seealso
#' \code{\link{CalculateDiffusionTF}}\cr
#' \code{vignette("calculate-transfer-functions")}
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#' 2018.
"diffusion.tf"

#' Time uncertainty transfer functions
#'
#' A dataset containing the calculated time uncertainty transfer functions for
#' the DML1, DML2 and WAIS records from Münch and Laepple (2018) based on
#' numerical simulations. The transfer function describes the effect time
#' uncertainty has on the spectrum of the average (in the time domain) of
#' a number of `N` proxy records.
#'
#' @format A named list with the 3 elements \code{dml1}, \code{dml2}, and
#'   \code{wais}, where each element is a spectral object:
#'   \describe{
#'   \item{\code{freq}:}{the frequency axis in units of \code{1 / yr},}
#'   \item{\code{spec}:}{the transfer function value at each frequency.}
#' }
#'
#' @source
#' The transfer functions were obtained using
#' \code{?CalculateTimeUncertaintyTF}; see also the respective package
#' vignette \code{vignette("calculate-transfer-functions")}.
#' @seealso
#' \code{\link{CalculateTimeUncertaintyTF}}\cr
#' \code{vignette("calculate-transfer-functions")}
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#' 2018.
"time.uncertainty.tf"

#' "proxysnr" spectral object
#'
#' @description{
#' A spectral object is a named list of the two elements `freq` and `spec`,
#' which are numeric vectors of the same length holding a frequency axis and the
#' corresponding power spectral densities.
#'
#' Optionally, the list can in addition include a vector with the degrees of
#' freedom of the spectrum (list element \code{dof}) and be of class "spec", but
#' "proxysnr" does not include any methods generic to that class as of now.
#' }
#' @name spec.object
NULL
