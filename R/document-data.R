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
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
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
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
"wais"

#' Trench oxygen isotope data
#'
#' A dataset containing the oxygen isotope data of 22 shallow firn
#' (\code{"trench"}) profiles from Kohnen Station, Dronning Maud Land,
#' Antarctica, as used in Münch and Laepple (2018).
#'
#' @format A data frame with 107 rows and 22 variables where each variable is the
#' depth profile of oxygen isotope data from one trench profile in per mil.
#' Vertical (depth) resolution is 3 cm.
#' @details
#' Data has been processed from the original source data as documented in the
#' package source (\code{data-raw/load.Trench.R}).
#' @source
#' Data citation:\cr
#' Münch, T., et al.: Stable water isotopes measured along two snow trenches
#' sampled at Kohnen Station, Dronning Maud Land, Antarctica in the 2014/15
#' field season, \url{https://doi.org/10.1594/PANGAEA.876639}, 2017.
#'
#' Literature citation:\cr
#' Münch, T., et al.: Constraints on post-depositional isotope modifications in
#' East Antarctic firn from analysing temporal changes of isotope profiles. The
#' Cryosphere, 11(5), 2175-2188, \url{https://doi.org/10.5194/tc-11-2175-2017},
#' 2017.\cr
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
"t15"

#' Climatological site parameters
#'
#' A dataset containing the climatological site parameters that are relevant for
#' the diffusion length calculation for all firn-core sites used in Münch and
#' Laepple (2018).
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{core:}{the name of the firn-core site}
#'   \item{temperature:}{10-m firn temperature (or annual mean air
#'     temperature), in degree Celsius}
#'   \item{acc.rate:}{accumulation rate, in kg/m^2/yr}
#'   \item{pressure:}{surface air pressure, in mbar}
#'   \item{height:}{surface elevation, in m}
#' }
#' @source
#' All values for firn cores in the rows 1-15 from Oerter et al. (2000) (Table
#' 1, 4), except: temperature for FB9807 adapted from B32 (1 km), temperature
#' for FB9815 adapted from FB9805 (28 km), temperature for FB9811 adapted from
#' FB9812 (19 km).
#'
#' For firn cores in the rows 16-20, accumulation rates and elevation data are
#' for core WDC2005A taken from WAIS Divide Project Members (2013) (p.440) and
#' for the other cores from Kaspari et al. (2004) (Table 1), respectively.
#' WDC2005A temperature is taken from WAIS Divide Project Members (2013)
#' (p.440), the other temperature values are based on ERA-Interim mean
#' temperature anomalies (Dee et al., 2011) with respect to WDC2005A.
#'
#' For all cores, the pressure values are calculated with the barometric formula
#' using the local elevation and temperature data.
#' @references
#' Dee, D. P., et al.: The ERA-Interim reanalysis: configuration and performance
#' of the data assimilation system, Q. J. R. Meteorol. Soc., 137(656), 553–597,
#' \url{https://doi.org/10.1002/qj.828}, 2011.
#' 
#' Kaspari, S., et al.: Climate variability in West Antarctica derived from
#' annual accumulation-rate records from ITASE firn/ice cores, Ann. Glaciol.,
#' 39(1), 585-594, \url{https://doi.org/10.3189/172756404781814447}, 2004.
#'
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
#'
#' Oerter, H., et al.: Accumulation rates in Dronning Maud Land, Antarctica, as
#' revealed by dielectric-profiling measurements of shallow firn cores,
#' Ann. Glaciol., 30(1), 27-34, 2000.
#'
#' WAIS Divide Project Members: Onset of deglacial warming in West Antarctica
#' driven by local orbital forcing, Nature, 500(7463), 440–444,
#' \url{https://doi.org/10.1038/nature12376}, 2013.
"clim.site.par"

#' Diffusion lengths
#'
#' A dataset containing calculated diffusion lengths for all firn-core sites
#' used in Münch and Laepple (2018). The diffusion lengths are provided as a
#' function of firn-core age at semi-annual resolution and in units of years.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{dml1:}{a data frame with 338 rows and 16 variables:}
#'     \describe{
#'       \item{Time:}{age of the firn core, in yr}
#'       \item{columns "B31"-"FB9817":}{15 columns with 338 rows with diffusion
#'       lengths corresponding to the firn-core ages in \code{Time}, in yr}
#'   }
#'   \item{dml2:}{a data frame with 1990 rows and 4 variables:}
#'     \describe{
#'       \item{Time:}{age of the firn core, in yr}
#'       \item{columns "B31"-"B33":}{3 columns with 1990 rows with diffusion
#'       lengths corresponding to the firn-core ages in \code{Time}, in yr}
#'   }
#'   \item{wais:}{a data frame with 402 rows and 6 variables:}
#'     \describe{
#'       \item{Time:}{age of the firn core, in yr}
#'       \item{columns "WDC2005A"-"ITASE-2000-5":}{5 columns with 402 rows with
#'       diffusion lengths corresponding to the firn-core ages in \code{Time},
#'       in yr}
#'   }
#' }
#' @details
#' Diffusion lengths as a function of firn-core depth are calculated from
#' Eq. (8) in Gkinis et al. (2014); firn diffusivity is modelled following
#' Eqs. (17-19) in Johnsen et al. (2000), firn density with the Herron–Langway
#' model (Herron and Langway, 1980; Eqs. 7 and 10). For the surface firn density
#' we assume a constant value of 340 kg/m^3 for all sites, for all other
#' relevant site parameters, see \code{?clim.site.par}. The depth-dependent
#' diffusion length values are converted into units of yr using the
#' Herron-Langway firn density with a constant accumulation rate, related to
#' the firn-core age and then linearly interpolated onto an equidistant time
#' axis at semi-annual resolution. All steps are documented in the package
#' source (\code{data-raw/calculate-diffusion-length.R}).
#' @references
#' Gkinis, V., et al.: Water isotope diffusion rates from the North-GRIP ice
#' core for the last 16,000 years – Glaciological and paleoclimatic
#' implications, Earth Planet. Sci. Lett., 405, 132–141,
#' \url{https://doi.org/10.1016/j.epsl.2014.08.022}, 2014.
#'
#' Herron, M. M. and Langway, Jr., C. C.: Firn Densification: An Empirical
#' Model, J. Glaciol., 25(95), 373–385,
#' \url{https://doi.org/10.3189/S0022143000015239}, 1980.
#'
#' Johnsen, S. J., et al.: Diffusion of stable isotopes in polar firn and ice:
#' the isotope effect in firn diffusion, in: Physics of Ice Core Records, edited
#' by Hondoh, T., vol. 159, pp. 121–140, Hokkaido University Press, Sapporo,
#' Japan, 2000.
#'
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
"diffusion.length"

#' Calculated diffusion transfer functions
#'
#' A dataset containing the calculated diffusion transfer functions for the
#' DML1, DML2 and WAIS data sets from Münch and Laepple (2018) based on
#' numerical simulations. The transfer function describes the average effect
#' diffusion has on the mean of the spectra from \code{N} spatially distributed
#' oxygen isotope records.
#'
#' @format A list with 3 elements:
#' \tabular{rll}{
#'   \code{dml1}: \tab a list with 2 numeric vectors of length 97: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#'   \code{dml2}: \tab a list with 2 numeric vectors of length 497: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#'   \code{wais}: \tab a list with 2 numeric vectors of length 100: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#' }
#' @source
#' The transfer functions were obtained using \code{?DiffusionTF} for the
#' site-specific diffusion lengths provided by \code{?diffusion.length}. Isotope
#' records were simulated at semiannual resolution and the transfer functions
#' interpolated in frequency space to annual resolution. See also the respective
#' package vignette:
#' \code{vignette(topic = "calculate-transfer-functions", package =
#'     "proxysnr")}.
#' @seealso
#' \code{\link{DiffusionTF}}\cr
#' \code{vignette(topic = "calculate-transfer-functions", package = "proxysnr")}
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
"diffusion.tf"

#' Calculated time uncertainty transfer functions
#'
#' A dataset containing the calculated time uncertainty transfer functions for
#' the DML1, DML2 and WAIS data sets from Münch and Laepple (2018) based on
#' numerical simulations. The transfer function describes the effect time
#' uncertainty has on the spectrum of the average (in the time domain) of
#' \code{N} proxy records.
#'
#' @format A list with 3 elements:
#' \tabular{rll}{
#'   \code{dml1}: \tab a list with 2 numeric vectors of length 97: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#'   \code{dml2}: \tab a list with 2 numeric vectors of length 497: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#'   \code{wais}: \tab a list with 2 numeric vectors of length 100: \tab  \cr
#'     \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
#'     \tab \code{spec}: transfer function value at each frequency \tab  \cr
#' }
#' @source
#' The transfer functions were obtained using \code{?TimeUncertaintyTF}. See
#' also the respective package vignette:
#' \code{vignette(topic = "calculate-transfer-functions", package =
#'     "proxysnr")}.
#' @seealso
#' \code{\link{DiffusionTF}}\cr
#' \code{vignette(topic = "calculate-transfer-functions", package = "proxysnr")}
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal to centennial scale isotope variations from Antarctic ice cores?
#' Clim. Past Discuss., \url{https://doi.org/10.5194/cp-2018-112}, in review,
#' 2018.
"time.uncertainty.tf"

