% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WrapSpectralResults.R
\name{WrapSpectralResults}
\alias{WrapSpectralResults}
\title{Wrap spectral results for several datasets}
\usage{
WrapSpectralResults(
  ...,
  diffusion = NULL,
  time.uncertainty = NULL,
  res = 1,
  df.log = 0.05,
  crit.diffusion = 0.5
)
}
\arguments{
\item{...}{a comma separated list of named proxy datasets to analyse.}

\item{diffusion}{a list the same length as the number of datasets with each
list element a spectral object (\code{?spec.object}) of a transfer function
to correct the corresponding dataset for the effect of diffusion-like
smoothing (\code{?SeparateSignalFromNoise} for more details on this, and
\code{?CalculateDiffusionTF} for calculating transfer functions
specifically for the firn diffusion process). Internally, the inverse of
the transfer function values are applied to correct for the smoothing
effect on the estimated signal and noise spectra (see Eq. 4 in Münch and
Laepple, 2018). If \code{NULL}, no correction is applied at all; if instead
you want to omit the correction only for some specific data set(s), set the
corresponding list element(s) to \code{NA}.}

\item{time.uncertainty}{as \code{diffusion} a list of transfer functions to
correct for the effect of time uncertainty
(\code{?CalculateTimeUncertaintyTF} for calculating transfer
functions in the case of layer-counted proxy chronologies).}

\item{res}{the sampling (e.g., temporal) resolution of the proxy data. Either
a single value if all datasets have the same resolution, or a vector with a
value for each dataset.}

\item{df.log}{a vector of Gaussian kernel widths in log space to smooth the
spectral estimates from each dataset. Either a single value to apply the
same smoothing to all datasets, or a vector with a value for each dataset.}

\item{crit.diffusion}{minimum transfer function value for the diffusion-like
smoothing process to constrain the corresponding correction. This
determines a cutoff frequency until which results are analysed to avoid
large uncertainties at the high-frequency end of the spectra; defaults to
0.5.}
}
\value{
A list of \code{n} lists, where \code{n} is the number of provided
  datasets and where each of these lists contains up to four elements:
  \describe{
  \item{\code{raw}:}{a list with four elements: three spectral objects (the
    raw signal, noise and corresponding SNR spectra), and a two-element
    vector (\code{f.cutoff}) with the index and value of the cutoff frequency
    from constraining the smoothing correction (see the \code{crit.diffusion}
    parameter).}
  \item{\code{corr.diff.only}:}{as item \code{raw} but with the spectra after
    correction for the effect of diffusion-like smoothing.}
  \item{\code{corr.t.unc.only}:}{as item \code{raw} but with the spectra
    after correction for the effect of time uncertainty.}
  \item{\code{corr.full}:}{as item \code{raw} but with the spectra after
    correction for both the effects of diffusion-like smoothing and time
    uncertainty.}
}
The number of the returned list elements for each dataset depends on
whether transfer functions for the corrections have been provided in
\code{diffusion} and \code{time.uncertainty} or not. Also, the element
\code{f.cutoff} is \code{NA} if diffusion-like smoothing has not been corrected
for.
}
\description{
This wrapper function is used to combine all main spectral results for the
signal, noise and signal-to-noise ratio as presented and discussed in Münch
and Laepple (2018), but it can also be used to combine the results for other
data sets.
}
\examples{
# Get main results of Münch and Laepple (2018)

results <- WrapSpectralResults(
               dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
               diffusion = diffusion.tf,
               time.uncertainty = time.uncertainty.tf,
               df.log = c(0.15, 0.15, 0.1))
}
\references{
Münch, T. and Laepple, T.: What climate signal is contained in
decadal- to centennial-scale isotope variations from Antarctic ice cores?
Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
  2018.
}
\seealso{
\code{\link{SeparateSignalFromNoise}},
  \code{\link{CalculateDiffusionTF}},
  \code{\link{CalculateTimeUncertaintyTF}}, \code{\link{spec.object}} for the
  definition of a \code{proxysnr} spectral object.
}
\author{
Thomas Münch
}
