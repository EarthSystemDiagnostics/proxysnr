% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PlotSNR}
\alias{PlotSNR}
\title{Plot proxy signal-to-noise ratios}
\usage{
PlotSNR(
  spec,
  f.cut = FALSE,
  names = NULL,
  col = 1:length(spec),
  xlim = c(500, 2),
  ylim = c(0.05, 5),
  xlab = "Time period (yr)",
  ylab = "Signal-to-Noise Ratio",
  xtm = NULL,
  ytm = NULL,
  xtl = NULL,
  ytl = NULL
)
}
\arguments{
\item{spec}{a (named) list of signal-to-noise ratio data sets: each dataset
itself should be list containing at least the spectral object
(\code{?spec.object}) \code{snr} providing signal-to-noise ratios as a
function of frequency. For Figure 3 in Münch and Laepple (2018) set
\code{spec} to the output from \code{\link{PublicationSNR}}.}

\item{f.cut}{Shall the spectra be cut at the cutoff frequency constrained
by the diffusion correction strength? Defaults to \code{FALSE}.}

\item{names}{an optional character vector of names of the proxy data
sets. If \code{NULL}, the names of \code{spec} are used or, if not present,
default names.}

\item{col}{a numeric or character vector of colors to use for the plotting
with length recycled to match \code{length(spec)}.}

\item{xlim}{the x limits (x1, x2) of the plot.}

\item{ylim}{the y limits (y1, y2) of the plot.}

\item{xlab}{x axis label.}

\item{ylab}{y axis label.}

\item{xtm}{x axis tick mark positions; if \code{NULL} computed by
\code{\link[graphics]{axis}}.}

\item{ytm}{y axis tick mark positions; default setting (\code{NULL}) uses
\code{c(0.05, 0.1, 0.5, 1, 5)}.}

\item{xtl}{x axis tick mark labels; if \code{NULL} determined automatically
from \code{xtm}, else it must be a vector of labels of the same length as
\code{xtm}.}

\item{ytl}{equivalent to \code{xtl} for the y axis tick mark labels.}
}
\description{
Plot proxy signal-to-noise ratios of several datasets as a function of
timescale (e.g., as in the firn core analysis of Münch and Laepple, 2018,
Fig. 3).
}
\examples{
# create toy data
n <- 100
spec <- list(
  data1 = list(snr = list(freq = seq(0.01, 0.5, length.out = n),
                          spec = seq(1, 0.1, length.out = n))),
  data2 = list(snr = list(freq = seq(0.005, 0.5, length.out = n),
                          spec = seq(5, 0.1, length.out = n)))
)

# plot SNR data
PlotSNR(spec)
}
\references{
Münch, T. and Laepple, T.: What climate signal is contained in
decadal- to centennial-scale isotope variations from Antarctic ice cores?
Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
  2018.
}
\seealso{
\code{\link{spec.object}} for the definition of a \code{proxysnr}
  spectral object.
}
\author{
Thomas Münch
}
