##
## Collection of plotting functions for the results in Münch and Laepple (2018)
##

#' Final DML and WAIS spectra
#'
#' Internal function to produce the final DML and WAIS spectra (signal, noise,
#' SNR) used for Figs. (3) and (4) in Münch and Laepple (2018). The DML spectra
#' are a combination of the results from the DML1 and DML2 data sets.
#'
#' @param spec spectral results for the DML and WAIS firn core data as output
#'   from \code{\link{WrapSpectralResults}}.
#' @param data character string naming the version of the spectral results in
#'   \code{spec} to use; one of "corr.full", "corr.diff.only",
#'   "corr.t.unc.only", and "raw"; defaults to "corr.full".
#' @param dml.knit.f frequency at which to combine the spectra from the DML1
#'   and DML2 data sets (defaults to 1/10 yr^(-1)); DML2 spectra are used for
#'   lower frequencies than \code{dml.knit.f}, DML1 for the higher frequencies.
#' @param df.log width of Gaussian kernel in logarithmic frequency units for
#'   additional smoothing for visual purposes; \code{NULL} suppresses smoothing.
#'
#' @return A list with the components \code{dml} and \code{wais}, where each is
#'   a list containing three spectral objects (`?spec.object`) and one numeric
#'   vector:
#'   \describe{
#'   \item{\code{signal}:}{the signal spectrum.}
#'   \item{\code{noise}:}{the noise spectrum.}
#'   \item{\code{snr}:}{the frequency-dependent signal-to-noise ratio.}
#'   \item{\code{f.cutoff}:}{a two-element vector with the index and value of
#'     the cutoff frequency from constraining the diffusion correction.}
#' }
#'
#' @author Thomas Münch
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#'   2018.
#'
#' @keywords internal
#'
PublicationSNR <- function(spec,
                           data = c("corr.full", "corr.diff.only",
                                    "corr.t.unc.only", "raw"),
                           dml.knit.f = 0.1, df.log = 0.125) {

  # Get correction version

  names.spec <- names(spec)
  nm <- match.arg(data,
                  c("corr.full", "corr.diff.only", "corr.t.unc.only", "raw"))

  for (i in 1 : length(spec)) {

    if (!all(utils::hasName(spec[[i]], nm))) {
      stop(paste0("No version `", nm, "` available for dataset `",
                  names.spec[i], "`."), call. = FALSE)
    }
  }

  spec.dml1 <- spec$dml1[[nm]]
  spec.dml2 <- spec$dml2[[nm]]
  spec.wais <- spec$wais[[nm]]


  # Combine DML1 and DML2 spectra

  idx.knit1 <- which(spec.dml1$signal$freq > dml.knit.f)
  idx.knit2 <- which(spec.dml2$signal$freq < dml.knit.f)

  dml.signal <- list()
  dml.signal$freq <- c(spec.dml2$signal$freq[idx.knit2],
                       spec.dml1$signal$freq[idx.knit1])
  dml.signal$spec <- c(spec.dml2$signal$spec[idx.knit2],
                       spec.dml1$signal$spec[idx.knit1])

  dml.noise <- dml.signal
  dml.noise$spec <- c(spec.dml2$noise$spec[idx.knit2],
                      spec.dml1$noise$spec[idx.knit1])

  dml.f.cutoff <- spec.dml1$f.cutoff
  if (length(dml.f.cutoff) == 2) {
    dml.f.cutoff[1] <- which(dml.signal$freq == spec.dml1$f.cutoff[2])
  }


  # Smooth DML & WAIS signal and noise and calculate SNR

  if (!is.null(df.log)) {

    dml.signal <- LogSmooth(dml.signal, df.log = df.log)
    dml.noise  <- LogSmooth(dml.noise, df.log = df.log)

    wais.signal <- LogSmooth(spec.wais$signal, df.log = df.log)
    wais.noise  <- LogSmooth(spec.wais$noise, df.log = df.log)
  }

  dml.snr <- list()
  dml.snr$freq <- dml.signal$freq
  dml.snr$spec <- dml.signal$spec / dml.noise$spec

  wais.snr <- list()
  wais.snr$freq <- wais.signal$freq
  wais.snr$spec <- wais.signal$spec / wais.noise$spec


  # Organize output

  res <- list()

  class(dml.signal) <- "spec"
  class(dml.noise)  <- "spec"
  class(dml.snr)    <- "spec"

  class(wais.signal) <- "spec"
  class(wais.noise)  <- "spec"
  class(wais.snr)    <- "spec"

  res$dml$signal   <- dml.signal
  res$dml$noise    <- dml.noise
  res$dml$snr      <- dml.snr
  res$dml$f.cutoff <- dml.f.cutoff

  res$wais$signal   <- wais.signal
  res$wais$noise    <- wais.noise
  res$wais$snr      <- wais.snr
  res$wais$f.cutoff <- spec.wais$f.cutoff

  attr(res$dml, "array.par") <- attr(spec.dml2, "array.par") # use DML2 parameter
  attr(res$wais, "array.par") <- attr(spec.wais, "array.par")

  return(res)

}

#' Figure 2 in Münch and Laepple (2018)
#'
#' Plot the signal and noise spectra from the DML and WAIS oxygen isotope data
#' sets (Figure 2 in Münch and Laepple, 2018).
#'
#' @param spec output from \code{\link{WrapSpectralResults}}.
#' @param f.cut Shall the spectra be cut at the cutoff frequency constrained
#'   by the diffusion correction strength? Defaults to \code{TRUE}.
#'
#' @author Thomas Münch
#' @seealso \code{\link{WrapSpectralResults}}
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#'   2018.
#'
#' @noRd
#'
muench_laepple_fig02 <- function(spec, f.cut = TRUE) {


  # --------------------------------------------------------------------------
  # Graphics settings
  
  ylabel <- expression("Power spectral density " * "(\u2030"^{2}%.%"yr)")
  ylim <- c(5e-2, 1e1)
  y.at <- c(0.05, 0.1, 0.5, 1, 5)
  removeLast <- 1

  op <- graphics::par(mar = c(0, 0, 0, 0), las = 1,
                      oma = c(5, 10, 2, 0.5), mfcol = c(2, 2),
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  on.exit(graphics::par(op))

  
  # --------------------------------------------------------------------------
  # Plot DML signal spectra

  LPlot(spec$dml1$raw$signal, type = "n", inverse = TRUE, axes = FALSE,
        xlim = c(500, 2), ylim = ylim, xlab = "", ylab = "")
  graphics::axis(2, at = y.at, labels = y.at)
  graphics::box()

  graphics::mtext(ylabel, side = 2, line = 4.5, las = 0,
                  cex = graphics::par()$cex.lab)
  graphics::mtext("DML", side = 2, line = 8, las = 0,
                  cex = graphics::par()$cex.lab)
  graphics::mtext("Signal", side = 3, line = 0.5, las = 0, adj = 0.99,
                  padj = 0.3, col = "dodgerblue4",
                  cex = graphics::par()$cex.lab)
  graphics::mtext("a", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)

  if (f.cut)
    removeLast <- length(
      spec$dml1$corr.full$f.cutoff[1] : length(spec$dml1$raw$signal$freq))
  
  LLines(spec$dml1$raw$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 1.5, lty = 3)
  LLines(spec$dml1$corr.t$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 1.5, lty = 5)
  LLines(spec$dml1$corr.full$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 3, lty = 1)

  if (f.cut)
    removeLast <- length(
      spec$dml2$corr.full$f.cutoff[1] : length(spec$dml2$raw$signal$freq))
  
  LLines(spec$dml2$raw$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 1.5, lty = 3)
  LLines(spec$dml2$corr.t$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 1.5, lty = 5)
  LLines(spec$dml2$corr.full$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 3, lty = 1)

  graphics::legend("bottomleft", c("DML1", "DML2"), lty = 1, lwd = 2,
                   col = c("dodgerblue4", "black"),  seg.len = 3,
                   bty = "n", cex = 1.25)
  graphics::legend("bottomleft",
                   c("Uncorrected signal", "Corrected for time uncertainty",
                     "Corrected for time uncertainty + diffusion"),
                   lwd = c(1.5, 1.5, 2), lty = c(3, 5, 1), col = "darkgrey",
                   inset = c(0.225, 0), seg.len = 3, bty = "n", cex = 1.25)

  # --------------------------------------------------------------------------
  # Plot WAIS signal spectra
  
  LPlot(spec$wais$raw$signal, type = "n", inverse = TRUE, axes = FALSE,
        xlim = c(500, 2), ylim = ylim, xlab = "", ylab = "")
  graphics::axis(1)
  graphics::axis(2, at = y.at, labels = y.at)
  graphics::box()

  graphics::mtext("Time period (yr)", side = 1, line = 3.5,
                  cex = graphics::par()$cex.lab)
  graphics::mtext(ylabel, side = 2, line = 4.5, las = 0,
                  cex = graphics::par()$cex.lab)
  graphics::mtext("WAIS", side = 2, line = 8, las = 0,
                  cex = graphics::par()$cex.lab)
  graphics::mtext("c", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)
  
  if (f.cut)
    removeLast <- length(
      spec$wais$corr.full$f.cutoff[1] : length(spec$wais$raw$signal$freq))

  LLines(spec$wais$raw$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 1.5, lty = 3)
  LLines(spec$wais$corr.t$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 1.5, lty = 5)
  LLines(spec$wais$corr.full$signal, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "dodgerblue4", lwd = 3, lty = 1)


  #---------------------------------------------------------------------------
  # Plot DML noise spectra
  
  LPlot(spec$dml1$raw$noise, type = "n", inverse = TRUE, axes = FALSE,
        xlim = c(500, 2), ylim = ylim, xlab = "", ylab = "")
  graphics::box()

  graphics::mtext("Noise", side = 3, line = 0.5, las = 0, adj = 0.015,
                  padj = 0.3, col = "firebrick4",
                  cex = graphics::par()$cex.lab)
  graphics::mtext("b", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)

  if (f.cut)
    removeLast <- length(
      spec$dml1$corr.full$f.cutoff[1] : length(spec$dml1$raw$signal$freq))
  
  LLines(spec$dml1$raw$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 1.5, lty = 3)
  LLines(spec$dml1$corr.t$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 1.5, lty = 5)
  LLines(spec$dml1$corr.full$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 3, lty = 1)

  if (f.cut)
    removeLast <- length(
      spec$dml2$corr.full$f.cutoff[1] : length(spec$dml2$raw$signal$freq))
  
  LLines(spec$dml2$raw$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 1.5, lty = 3)
  LLines(spec$dml2$corr.t$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 1.5, lty = 5)
  LLines(spec$dml2$corr.full$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "black", lwd = 3, lty = 1)

  graphics::legend("bottomleft", c("DML1", "DML2"), lty = 1, lwd = 2,
                   col = c("firebrick4", "black"),  seg.len = 3,
                   bty = "n", cex = 1.25)
  graphics::legend("bottomleft",
                   c("Uncorrected noise", "Corrected for time uncertainty",
                     "Corrected for time uncertainty + diffusion"),
                   lwd = c(1.5, 1.5, 2), lty = c(3, 5, 1), col = "darkgrey",
                   inset = c(0.225, 0), seg.len = 3, bty = "n", cex = 1.25)

  #---------------------------------------------------------------------------
  # Plot WAIS noise spectra
  
  LPlot(spec$wais$raw$noise, type = "n", inverse = TRUE, axes = FALSE,
        xlim = c(500, 2), ylim = ylim, xlab = "", ylab = "")
  graphics::axis(1)
  graphics::box()

  graphics::mtext("Time period (yr)", side = 1, line = 3.5,
                  cex = graphics::par()$cex.lab)
  graphics::mtext("d", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)

  if (f.cut)
    removeLast <- length(
      spec$wais$corr.full$f.cutoff[1] : length(spec$wais$raw$signal$freq))

  LLines(spec$wais$raw$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 1.5, lty = 3)
  LLines(spec$wais$corr.t$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 1.5, lty = 5)
  LLines(spec$wais$corr.full$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 3, lty = 1)

}

#' Figure 5 in Münch and Laepple (2018)
#'
#' Plot comparison of the DML and T15 trench noise spectra (Figure 5 in Münch
#' and Laepple, 2018).
#'
#' @param SNR output from \code{\link{PublicationSNR}}.
#' @param TNS T15 trench noise spectra; per default loaded from internal data.
#' @param f.cut Shall the spectra be cut at the cutoff frequency constrained
#'   by the diffusion correction strength? Defaults to \code{TRUE}.
#'
#' @author Thomas Münch
#' @seealso \code{\link{PublicationSNR}}
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#'   2018.
#'
#' @noRd
#'
muench_laepple_fig05 <- function(SNR, TNS = t15.noise, f.cut = TRUE) {

  # Graphics settings

  xlab = "Time period (yr)"
  xlim = c(50, 0.5)
  xtm <- c(50, 20, 10, 5, 2, 1, 0.5)
  ylab <- expression("Noise PSD " * "(\u2030"^{2}%.%"yr)")
  ylim <- c(5e-2, 1e1)
  ytm <- c(0.05, 0.1, 0.5, 1, 5, 10)

  removeLast <- 0
  if (f.cut)
    removeLast <-
      length(SNR$dml$f.cutoff[1] : length(SNR$dml$noise$freq))

  op <- graphics::par(mar = c(5, 6.5, 0.5, 0.5), las = 1,
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
  on.exit(graphics::par(op))

  # Plot final DML noise spectrum

  LPlot(SNR$dml$noise, type = "n", inverse = TRUE, axes = FALSE,
        xlim = xlim, ylim = ylim, xlab = "", ylab = "")

  LLines(SNR$dml$noise, conf = FALSE, inverse = TRUE,
         removeFirst = 1, removeLast = removeLast,
         col = "firebrick4", lwd = 3, lty = 1)

  # Plot trench noise spectrum
  
  i.remove <- c(1, 2, length(TNS$lower$freq))

  # shaded spectral range according to upper/lower accumulation rate
  graphics::polygon(
              c(
                1 / TNS$lower$freq[-i.remove],
                rev(1 / TNS$upper$freq[-i.remove])),
              c(
                TNS$lower$spec[-i.remove],
                rev(TNS$upper$spec[-i.remove])
              ),
              col = grDevices::adjustcolor("dodgerblue4", 0.2), border = NA)

  # trench noise spectrum for mean accumulation rate
  LLines(TNS$mean, conf = FALSE, inverse = TRUE,
         removeFirst = 2, removeLast = 0,
         col = "dodgerblue4", lwd = 3, lty = 1)

  # Axis and legends settings

  graphics::axis(1, at = xtm, labels = xtm)
  graphics::axis(2, at = ytm, labels = ytm)

  graphics::mtext(xlab, side = 1, line = 3.5,
                  cex = graphics::par()$cex.lab)
  graphics::mtext(ylab, side = 2, line = 4.5, las = 0,
                  cex = graphics::par()$cex.lab)

  graphics::legend("bottomleft", c("Array scale (DML data set)",
                                   "Local scale (trench data set)"),
                   col = c("firebrick4", "dodgerblue4"),
                   seg.len = 3, lty = 1, lwd = 2, bty = "n")

  graphics::par(op)

}

