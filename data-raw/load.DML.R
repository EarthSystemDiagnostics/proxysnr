## load.DML.R: Download and process DML oxygen isotope data
##
## Author: Thomas Münch (thomas.muench@awi.de), Alfred-Wegener-Institut, 2018
##


# 'pangaear' package for easy download from PANGAEA repository
library(pangaear)

# PANGAEA doi's of the DML data sets
doi <- c("10.1594/PANGAEA.104889",
         "10.1594/PANGAEA.104890",
         "10.1594/PANGAEA.104880",
         "10.1594/PANGAEA.104891",
         "10.1594/PANGAEA.104879",
         "10.1594/PANGAEA.104892",
         "10.1594/PANGAEA.104893",
         "10.1594/PANGAEA.104878",
         "10.1594/PANGAEA.104887",
         "10.1594/PANGAEA.104886",
         "10.1594/PANGAEA.104885",
         "10.1594/PANGAEA.104884",
         "10.1594/PANGAEA.104882",
         "10.1594/PANGAEA.104881",
         "10.1594/PANGAEA.104888")

# Names of the respective data sets
names <- c(paste("fb",
                 c("04", "05", "07", "08", "09", "10",
                 "11", "13", "14", "15", "16", "17"),
                 sep = ""),
           "b31", "b32", "b33")

# Download and process entire DML data
dml1 <- load.DML(doi = doi, names = names,
                 setStart = 1994,
                 setEnd = c(rep(1801, 12), rep(1000, 3)),
                 fillNA = TRUE,
                 clearCache = TRUE, verbose = TRUE)

# Split into the two DML data sets used in the paper
dml2 <- structure(dml1[13 : 15], info = attr(dml1, "info"))
dml1[13 : 15] <- lapply(dml1[13 : 15], function(x){x[1 : length(dml1[[1]])]})

dml <- list(dml1 = dml1, dml2 = dml2)


#-------------------------------------------------------------------------------
##' load.DML: Download and process DML data
##'
##' Download DML isotope data from the PANGAEA repository and process them for
##' the spectral analyses.
##' @param doi vector of PANGAEA doi's ("10.1594/PANGAEA.xxxxxx")
##' @param names optional character vector of the names of the data sets
##' @param setStart vector of start dates for each data set; if \code{NULL} the
##' original start date is used. If only one date is given, the length of
##' \code{setStart} is recycled to match the number of data sets.
##' @param setEnd vector of end dates for each data set; if \code{NULL} the
##' original end date is used. If only one date is given, the length of
##' \code{setEnd} is recycled to match the number of data sets.
##' @param fillNA Shall missing values be linearly interpolated? Defaults to
##' \code{FALSE}.
##' @param clearCache Shall the local pangaea data cache (the downloaded data
##' text files) be deleted? Defaults to \code{FALSE}.
##' @param verbose if \code{TRUE}, messages on the download and processing are
##' printed.
##' @return a list with the oxygen isotope data for the data sets.
##' @author Thomas Münch
load.DML <- function(doi,
                     names = NULL,
                     setStart = NULL,
                     setEnd = NULL,
                     fillNA = FALSE,
                     clearCache = FALSE,
                     verbose = FALSE) {

    if (length(setStart) == 1) {
        setStart <- rep(setStart, length(doi))
    }
    if (length(setEnd) == 1) {
        setEnd <- rep(setEnd, length(doi))
    }

    # Download data from PANGAEA
    data <- list()
    for (i in 1 : length(doi)) {

        tmp <- pangaear::pg_data(doi = doi[i], mssgs = verbose)
        
        # Keep age (yr AD) and d18O values only
        data[[i]] <- as.data.frame(
            tmp[[1]]$data)[, c("Age [a AD]", "d18O H2O [per mil SMOW]")]

        # Adjust time frame of record
        i.start <- ifelse(is.null(setStart[i]),
                          1, match(setStart[i], data[[i]][, 1]))
        i.end <- ifelse(is.null(setEnd[i]),
                        nrow(data[[i]]), match(setEnd[i], data[[i]][, 1]))

        data[[i]] <- data[[i]][i.start : i.end, ]

    }
    
    # Set names for the data set if provided
    if (!is.null(names)) {
        names(data) <- names
    }

    # Interpolate missing values if desired
    if (fillNA) {

        # Print number of interpolated vs. total number of data points
        if (verbose) {
            
            na.amount <- sapply(data, function(df) {
                sum(is.na(df[, 2]))
            })

            print(sprintf("Total # data points: %2.0f",
                          length(unlist(data)) / 2))
            print(sprintf("Data points interpolated: %2.0f.",
                          sum(na.amount)))

        }

        data <- lapply(data, function(df) {
            df[, 2] <- approx(df[, 1], df[, 2], df[, 1])$y
            return(df)
        })
        
    }

    # Keep only d18O values
    data <- lapply(data, function(df) {df[, 2]})
    
    # Clear the local pangaea data cache if desired
    if (clearCache) {
        pangaear::pg_cache_clear(prompt = FALSE)
    }

    # Provide meta info as attributes
    attr(data, "info") <-
        sprintf("Downloaded and processed on: %s.", Sys.time())

    return(data)

}

