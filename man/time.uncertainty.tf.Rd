% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/document-data.R
\docType{data}
\name{time.uncertainty.tf}
\alias{time.uncertainty.tf}
\title{Calculated time uncertainty transfer functions}
\format{
A list with 3 elements:
\tabular{rll}{
  \code{dml1}: \tab a list with 2 numeric vectors of length 97: \tab  \cr
    \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
    \tab \code{spec}: transfer function value at each frequency \tab  \cr
  \code{dml2}: \tab a list with 2 numeric vectors of length 497: \tab  \cr
    \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
    \tab \code{spec}: transfer function value at each frequency \tab  \cr
  \code{wais}: \tab a list with 2 numeric vectors of length 100: \tab  \cr
    \tab \code{freq}: frequency axis, in yr^{-1} \tab  \cr
    \tab \code{spec}: transfer function value at each frequency \tab  \cr
}
}
\source{
The transfer functions were obtained using \code{?TimeUncertaintyTF}. See
also the respective package vignette:
\code{vignette(topic = "calculate-transfer-functions", package =
    "proxysnr")}.
}
\usage{
time.uncertainty.tf
}
\description{
A dataset containing the calculated time uncertainty transfer functions for
the DML1, DML2 and WAIS data sets from Münch and Laepple (2018) based on
numerical simulations. The transfer function describes the effect time
uncertainty has on the spectrum of the average (in the time domain) of
\code{N} proxy records.
}
\references{
Münch, T. and Laepple, T.: What climate signal is contained in
decadal- to centennial-scale isotope variations from Antarctic ice cores?
Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
2018.
}
\seealso{
\code{\link{DiffusionTF}}\cr
\code{vignette(topic = "calculate-transfer-functions", package = "proxysnr")}
}
\keyword{datasets}
