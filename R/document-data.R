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

