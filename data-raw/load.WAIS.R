## load.WAIS.R: Download and process WAIS oxygen isotope data
##
## Author: Thomas Münch (thomas.muench@awi.de), Alfred-Wegener-Institut, 2018
##

library(usethis)
library(Hmisc)

# URL of the entire WAIS data sets
addr <- "http://www.usap-dc.org/dataset/nsidc/nsidc0536_steig/" 
url <- paste(addr, "Steig_2013_IsotopeAnnualAverages.txt", sep = "")

# Names of cores to select
cores <- c("WDC05A", "ITASE_1999_1",
           "ITASE_2000_1", "ITASE_2000_4", "ITASE_2000_5")

# Download and process WAIS data
wais <- load.WAIS(url = url,
                  cores = cores,
                  cores.seasonal2annual = cores[4 : 5],
                  fillNA = TRUE, verbose = TRUE)

usethis::use_data(wais, overwrite = TRUE)


#-------------------------------------------------------------------------------
# utility functions

#' load.WAIS: Download and process WAIS data
#'
#' Download WAIS isotope data from repository and process them for the spectral
#' analyses.
#'
#' @param url URL of the UASAP-DC repository data file with the WAIS annual
#'   averages data.
#' @param cores name of the WAIS cores to extract from the data file.
#' @param cores.seasonal2annual name of the cores for which annual data has to
#'   be obtained from high-resolution seasonal data (will be downloaded
#'   additionally).
#' @param setStart set common start date (year) for all cores; if `NULL` the
#'   original start date is used.
#' @param setEnd set common end date (year) for all cores; if `NULL` the
#'   original end date is used.
#' @param fillNA Shall missing values be linearly interpolated? Defaults to
#'   `FALSE`.
#' @param verbose if `TRUE`, messages on the processing are printed.
#'
#' @return a `data.frame` with the requested oxygen isotope data.
#'
#' @author Thomas Münch
#'
load.WAIS <- function(url, cores,
                      cores.seasonal2annual = NULL,
                      setStart = 2000, setEnd = 1800,
                      fillNA = FALSE,
                      verbose = FALSE) {

  seasonal2annual <- function(path, file, start, end) {

    url <- file.path(path, file)

    r <- readLines(url)
    i <- grep("% Depth_top", r)

    data <- read.table(url, skip = i - 1, header = TRUE,
                       sep = "\t", na.strings = "999999")

    # Mean age of data
    data$mean.age <- rowMeans(data[, c("Age_top", "Age_bottom")])

    # Bin-average to annual resolution
    d18O.annual <-
      rev(AvgToBin(rev(data[, "mean.age"]), rev(data[, "d18O..per.mil."]),
                   breaks = seq(end, start + 1, 1))$avg)

    return(d18O.annual)

  }

  # Header info has to be read explicitly for this file
  header <- read.table(url, skip = 26, nrow = 1, sep = "\t",
                       stringsAsFactors = FALSE)
  # Read actual data
  raw <- read.table(url, skip = 28, sep = "\t", na.strings = "999999")
  # Add proper header
  colnames(raw) <- c("Year", unlist(header)[-1])

  # Select cores and set time frame
  i.c <- match(cores, names(raw))
  i.t <- match(setStart : setEnd, raw$Year)
  data <- raw[i.t, c(1, i.c)]

  # Bin-average high-resolution data to annual resolution if needed
  if (!is.null(cores.seasonal2annual)) {

    path <- dirname(url)

    for (cr in cores.seasonal2annual) {

      data[, cr] <- seasonal2annual(path, paste(cr, "txt", sep = "."),
                                    setStart, setEnd)

    }
  }

  # Interpolate/Extrapolate missing values if desired
  if (fillNA) {

    # Print number of interpolated vs. total number of data points
    if (verbose) {

      na.amount <- sum(is.na(unlist(data[, -1])))

      print(sprintf("Total # data points: %2.0f",
                    length(unlist(data[, -1]))))
      print(sprintf("Data points interpolated: %2.0f.",
                    na.amount))

    }

    # Linearly interpolate missing values inside time series
    for (i in 1 : ncol(data)) {
      data[, i] <- approx(data$Year, data[, i], data$Year)$y
    }

    # Linearly extrapolate missing values at upper end of cores
    na.col <- which(is.na(data[1, ]))
    for (i in na.col) {
      j <- which(!is.na(data[, i][1 : 25]))
      data[1 : 25, i] <- Hmisc::approxExtrap(
                                  data$Year[j], data[j, i], data$Year[1 : 25])$y
    }
  }

  # Keep only d18O values
  data$Year <- NULL

  # Provide meta info as attributes
  attr(data, "info") <-
    sprintf("Downloaded and processed on: %s.", Sys.time())

  return(data)

}

#' Bin averaging
#'
#' Average a vector into bins.
#'
#' This function averages the vector `y` into bins according to the positon
#' of `x` within the breaks; you can either specify a desired number `N` of
#' breaks or directly specify the N + 1 break positions.
#'
#' @param x vector of values on which the data in `y` is tabulated; e.g. depth
#'   or time points.
#' @param y vector of observation values to be averaged into bins. Must have the
#'   same length as `x`.
#' @param N desired number of breaks (ignored if `breaks` are supplied
#'   directly).
#' @param breaks vector of break point positions to define the averagig bins; if
#'   omitted, break point positions are calculated from the range of `x` and the
#'   desired number of breaks given by `N`.
#' @param right logical; indicate whether the bin intervals should be closed on
#'   the right and open on the left (`TRUE`, the default), or vice versa
#'   (`FALSE`).
#' @param bFill logical; if `TRUE`, fill empty bins using linear interpolation
#'   from the neighbours to the center of the bin.
#'
#' @return a list with four elements:
#' \describe{
#' \item{`breaks`:}{numeric vector of the used break point positions.}
#' \item{`centers`:}{numeric vector with the positions of the bin centers.}
#' \item{`avg`:}{numeric vector with the bin-averaged values.}
#' \item{`nobs`:}{numeric vector with the number of observations
#'   contributing to each bin average.}
#' }
#'
#' @author Thomas Laepple
#'
AvgToBin <- function(x, y, N = 2, breaks = pretty(x, N),
                     right = TRUE, bFill = FALSE) {

  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }

  nBins <- length(breaks) - 1

  centers <- (breaks[1 : nBins] + breaks[2 : (nBins + 1)]) / 2
  nObs <- avg <- rep(NA, nBins)

  for (i in 1 : nBins) {

    if (right) {
      selection <- y[which((x > breaks[i]) & (x <= breaks[i + 1]))]
    } else {
      selection <- y[which((x >= breaks[i]) & (x < breaks[i + 1]))]
    }

    avg[i]  <- mean(na.omit(selection))
    nObs[i] <- sum(!is.na(selection))

  }

  if ((sum(missing <- is.na(avg)) > 0) & (bFill)) {

    avg[missing] <- (approx(x, y, centers)$y)[missing]

  }

  list(breaks = breaks, centers = centers, avg = avg, nobs = nObs)

}
