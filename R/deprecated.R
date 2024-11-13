#' Deprecated functions
#'
#' @description
#' The following function names are deprecated, since they do not adhere to the
#' paradigm of using verb phrases, or are not desriptive enough, and are
#' replaced by the respective new names mentioned thereafter:
#'
#' \itemize{
#'   \item \code{ArraySpectra} -> \code{\link{ObtainArraySpectra}}
#'   \item \code{SeparateSpectra} -> \code{\link{SeparateSignalFromNoise}}
#'   \item \code{StackCorrelation} -> \code{\link{ObtainStackCorrelation}}
#'   \item \code{IntegratedSNR} -> \code{\link{GetIntegratedSNR}}
#'   \item \code{DiffusionTF} -> \code{\link{CalculateDiffusionTF}}
#'   \item \code{TimeUncertaintyTF} -> \code{\link{CalculateTimeUncertaintyTF}}
#' }
#'
#' However, calls that use the old name will still work and pipe the input to
#' the respective function of the new name, but a warning will be issued
#' informing about the new function name.
#'
#' @name proxysnr-deprecated
NULL

#' @usage NULL
#' @export
ArraySpectra <- function(cores, res = 1, neff = length(cores),
                         df.log = NULL, ...) {

  .Deprecated("ObtainArraySpectra", "proxysnr")
  ObtainArraySpectra(cores, res, neff, df.log, ...)

}

#' @usage NULL
#' @export
SeparateSpectra <- function(spectra, neff = spectra$N,
                            diffusion, time.uncertainty) {

  if (missing(diffusion)) diffusion <- NULL
  if (missing(time.uncertainty)) time.uncertainty <- NULL

  .Deprecated("SeparateSignalFromNoise", "proxysnr")
  SeparateSignalFromNoise(spectra, neff, diffusion, time.uncertainty)

}

#' @usage NULL
#' @export
StackCorrelation <- function(input, N = 1, f1 = 2, f2 = "max",
                             freq.cut.lower = NULL, freq.cut.upper = NULL) {

  limits <- c(freq.cut.lower, freq.cut.upper)
    # <- does not work when only one of the two is specified, but best I can do

  .Deprecated("ObtainStackCorrelation", "proxysnr")
  ObtainStackCorrelation(input, N, f1, f2, limits)

}

#' @usage NULL
#' @export
IntegratedSNR <- function(input, N = 1, f1 = 2, f2 = "max",
                          freq.cut.lower = NULL, freq.cut.upper = NULL) {

  limits <- c(freq.cut.lower, freq.cut.upper)
    # <- does not work when only one of the two is specified, but best I can do

  .Deprecated("GetIntegratedSNR", "proxysnr")
  GetIntegratedSNR(input, N, f1, f2, limits)

}

#' @usage NULL
#' @export
DiffusionTF <- function(nt, nc, ns, sigma, res = 1, window = NULL,
                        coherent = FALSE, ...) {

  .Deprecated("CalculateDiffusionTF", "proxysnr")
  CalculateDiffusionTF(nt, nc, ns, sigma, res, window, coherent,
                       verbose.output = TRUE, ...)

}

#' @usage NULL
#' @export
TimeUncertaintyTF <- function(t = 100 : 1, acp = c(t[1], NA),
                              nt = length(t), nc = 1, ns = 100,
                              model = "poisson", rate = 0.05, resize = 1,
                              surrogate.fun = stats::rnorm, fun.par = NULL,
                              pad = TRUE, ...) {

  .Deprecated("CalculateDiffusionTF", "proxysnr")
  CalculateTimeUncertaintyTF(t, acp, nt, nc, ns, model, rate, resize,
                             surrogate.fun, fun.par, pad, verbose.output = TRUE,
                             ...)

}
