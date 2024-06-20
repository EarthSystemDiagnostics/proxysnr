##
## Collection of general plotting functions
##

# ------------------------------------------------------------------------------
# unexported utility functions

#' Draw error shading
#'
#' A wrapper function for the \code{polygon} function to draw error shadings
#' (or confidence intervals) on a line plot.
#'
#' @param x numeric vector of x values of the error band.
#' @param y1 numeric vector for the upper bound of the error band; must be of
#'   the same length as \code{x}.
#' @param y2 numeric vector for the lower bound of the error band; must be of
#'   the same length as \code{x}.
#' @param col colour of the error band.
#' @param alpha opacity factor for \code{col} within [0,1].
#' @param ... additional parameters which are passed to \code{polygon}.
#'
#' @author Thomas Münch
#' @noRd
#'
Polyplot <- function(x, y1, y2, col = "black", alpha = 0.2, ...) {

  inp <- list(x, y1, y2)
  if (stats::var(sapply(inp, length)) != 0)
    stop("All input vectors must be of the same length.")
  if (any(sapply(inp, function(x){any(is.na(x))})))
    warning("Polyplot: Missing values as input.", call. = FALSE)

  col <- grDevices::adjustcolor(col = col, alpha = alpha)

  graphics::polygon(c(x, rev(x)), c(y1, rev(y2)),
                    col = col, border = NA, ...)
  
}

#' Log-log spectral plot
#'
#' This function plots a spectrum on a double-logarithmic scale and optionally
#' adds a transparent confidence interval.
#'
#' @param x an object of class \code{"spec"}.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#'   interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} the x-axis is displayed in units of period
#'   (inverse frequency), increasing to the left. Defaults to \code{FALSE}.
#' @param bNoPlot if \code{TRUE} only produce the plot frame (\code{type = "n"}
#'   behaviour of function \code{\link{plot}}). Defaults to \code{FALSE}.
#' @param axes if \code{FALSE} the plotting of the x and y axes is
#'   suppressed. Defaults to \code{TRUE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#'   interval. Defaults to \code{0.2}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side.
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param xlab character string for labelling the x-axis.
#' @param ylab character string for labelling the y-axis.
#' @param xlim range of x-axis values; if \code{NULL} (the default) it is
#'   calculated internally and automatically reversed for \code{bPeriod = TRUE}.
#' @param ylim range of y-axis values; if \code{NULL} (the default) it is
#'   calculated internally.
#' @param ... further graphical parameters passed to \code{plot}.
#'
#' @author Thomas Laepple
#' @noRd
#'
LPlot <- function(x, conf = TRUE, bPeriod = FALSE, bNoPlot = FALSE, axes = TRUE,
                  col = "black", alpha = 0.2, removeFirst = 0, removeLast = 0,
                  xlab = "f", ylab = "PSD", xlim = NULL, ylim = NULL, ...) {

  if (bPeriod) {
    x$freq <- 1 / x$freq
    if (is.null(xlim)) xlim <- rev(range(x$freq))
    if (xlab == "f") xlab <- "period"
  }

  index <- (removeFirst + 1) : (length(x$freq) - removeLast)
  
  x$freq  <- x$freq[index]
  x$spec  <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]

  graphics::plot(x$freq, x$spec, type = "n", log = "xy",
                 xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                 axes = axes, ...)

  lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
  if (conf & lim & !bNoPlot) {
    Polyplot(x = x$freq, y1 = x$lim.1, y2 = x$lim.2,
             col = col, alpha = alpha)
  }

  if (!bNoPlot) graphics::lines(x$freq, x$spec, col = col, ...)
  
}

#' Add spectrum to existing log-log spectral plot
#'
#' This function adds a spectrum to an existing double-logarithmic plot and
#' optionally adds a transparent confidence interval.
#' 
#' @param x an object of class \code{"spec"}.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#'   interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} treat the x-axis values in units of period
#'   (inverse frequency). Defaults to \code{FALSE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#'   interval. Defaults to \code{0.2}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side. 
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param ... further graphical parameters passed to \code{lines}.
#'
#' @author Thomas Laepple
#' @noRd
#'
LLines<-function(x, conf = TRUE, bPeriod = FALSE, col = "black", alpha = 0.2,
                 removeFirst = 0, removeLast = 0, ...) {

  if (bPeriod) x$freq <- 1 / x$freq
  
  index <- (removeFirst + 1) : (length(x$freq) - removeLast)

  x$freq  <- x$freq[index]
  x$spec  <- x$spec[index]
  x$lim.1 <- x$lim.1[index]
  x$lim.2 <- x$lim.2[index]

  lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
  if (conf & lim) {
    Polyplot(x = x$freq, y1 = x$lim.1, y2 = x$lim.2,
             col = col, alpha = alpha)
  }
  
  graphics::lines(x$freq, x$spec, col = col, ...)
  
}

# ------------------------------------------------------------------------------
# generic functions to plot results from applying proxysnr method

#' Plot proxy array spectra
#'
#' Plot the spectral estimates from a spatial proxy core array, such as in
#' Fig. 1 in Münch and Laepple, 2018.
#'
#' @param spec output from \code{\link{ObtainArraySpectra}}.
#' @param f.cutoff optional frequency value to draw a vertical line symbolising
#'   a cutoff frequency.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim the y limits (y1, y2) of the plot.
#' @param col a three-element vector with the colours for the individual
#'   spectra (\code{col[1]}), the mean (\code{col[2]}) and the stack spectrum
#'   (\code{col[3]}).
#' @param col.sn a two-element vector with the colours for the approximate
#'   signal (\code{col[1]}) and noise shading (\code{col[2]}).
#' @param alpha.sn opacity factor for the colours in \code{col.sn} within
#'   [0,1].
#' @param plt.ann if \code{"default"} use axis annotation as in Fig. 1 of Münch
#'   and Laepple (2018). Since no other fixed annotation scheme is implemented,
#'   setting \code{plt.ann} to a different value will result in an error.
#' @param xlab if not \code{NULL} use this specific x axis label to override
#'   the default setting.
#' @param ylab if not \code{NULL} use this specific y axis label to override
#'   the default setting.
#' @param xtm if not \code{NULL} use this vector of specific x axis tick mark
#'   positions to override the default setting.
#' @param ytm if not \code{NULL} use this vector of specific y axis tick mark
#'   positions to override the default setting.
#' @param xtl x axis tick mark labels. If \code{xtm} is not \code{NULL},
#'   \code{xtl} is used for the \code{labels} parameter of the
#'   \code{\link{axis}} function: if \code{xtl} is set to \code{NULL}, the
#'   labels are determined automatically, else set it explicitly by specifying a
#'   vector of labels of the same length as \code{xtm}; see also \code{?axis}.
#' @param ytl equivalent to \code{xtl} for the y axis tick mark labels.
#'
#' @author Thomas Münch
#' @seealso \code{\link{ObtainArraySpectra}}
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @examples
#'
#' # Plot Figure 1 in Münch and Laepple (2018) (DML1 oxygen isotope data set):
#' PlotArraySpectra(ObtainArraySpectra(dml$dml1, df.log = 0.12))
#'
#' @export
#'
PlotArraySpectra <- function(spec, f.cutoff = NA,
                             xlim = c(100, 2), ylim = c(0.005, 50),
                             col = c("darkgrey", "black", "burlywood4"),
                             col.sn = c("dodgerblue4", "firebrick4"),
                             alpha.sn = 0.2,
                             plt.ann = "default",
                             xlab = NULL, ylab = NULL,
                             xtm = NULL, ytm = NULL,
                             xtl = NULL, ytl = NULL) {

  # Gather input
  
  N <- length(spec$single)    
  psd1 <- spec$mean
  psd2 <- spec$stack
  psd3 <- psd1
  psd3$spec <- psd3$spec / N

  # Axis settings

  if (plt.ann == "default") {
    set.xlab <- "Time period (yr)"
    set.ylab <- expression("Power spectral density " * "(\u2030"^{2}%.%"yr)")
    set.xtm <- NULL
    set.xtl <- NULL
    set.ytm <- 10^(seq(log10(ylim[1]), log10(ylim[2]), by = 1))
    set.ytl <- format(set.ytm,
                      scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  } else {
    stop("Invalid axis annotation setting.")
  }
  if (!is.null(xlab)) set.xlab = xlab
  if (!is.null(ylab)) set.ylab = ylab
  if (!is.null(xtm)) {set.xtm = xtm; set.xtl = xtl}
  if (!is.null(ytm)) {set.ytm = ytm; set.ytl = ytl}
  
  # Plot parameters

  op <- graphics::par(mar = c(5, 6.5, 0.5, 0.5), las = 1,
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
  on.exit(graphics::par(op))

  # Plot frame
  
  LPlot(psd1, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
        xlim = xlim, ylim = ylim, xlab = "", ylab = "")

  # Shadings
  
  i.remove <- c(1, length(psd1$freq))

  Polyplot(1 / psd1$freq[-i.remove],
           psd1$spec[-i.remove], psd2$spec[-i.remove],
           col = col.sn[2], alpha = alpha.sn)
  Polyplot(1 / psd1$freq[-i.remove],
           psd2$spec[-i.remove], psd3$spec[-i.remove],
           col = col.sn[1], alpha = alpha.sn)

  # Frequency cutoff line

  graphics::lines(x = rep(1 / f.cutoff, 2), y = c(ylim[1]/10, ylim[2]),
                  lty = 2, col = "darkgrey")

  # Individual spectra

  for (i in 1 : N) {
    
    LLines(spec$single[[i]], conf = FALSE, bPeriod = TRUE,
           removeFirst = 1, removeLast = 1,
           col = col[1])
  }

  # Main spectra

  LLines(psd1, conf = FALSE, bPeriod = TRUE,
         removeFirst = 1, removeLast = 1,
         col = col[2], lwd = 3)
  LLines(psd2, conf = FALSE, bPeriod = TRUE,
         removeFirst = 1, removeLast = 1,
         col = col[3], lwd = 3)
  LLines(psd3, conf = FALSE, bPeriod = TRUE,
         removeFirst = 1, removeLast = 1,
         col = col[2], lwd = 1.5, lty = 5)

  # Axis and legend settings

  graphics::axis(1, at = set.xtm, labels = set.xtl)
  graphics::axis(2, at = set.ytm, labels = set.ytl)

  graphics::mtext(set.xlab, side = 1, line = 3.5,
                  cex = graphics::par()$cex.lab)
  graphics::mtext(set.ylab, side = 2, line = 4.5,
                  cex = graphics::par()$cex.lab, las = 0)

  graphics::legend("bottomleft",
                   c("Individual spectra", "Mean spectrum",
                     "Spectrum of stacked record", "Mean spectrum scaled by 1/n"),
                   col = c(col[1 : 3], col[2]),
                   lty = c(1, 1, 1, 5),
                   lwd = c(1, 2, 2, 1), seg.len = 2.5, bty = "n")

}

#' Plot proxy signal-to-noise ratios
#'
#' Plot the timescale dependence of proxy signal-to-noise ratios, such as in
#' Fig. 3 in Münch and Laepple, 2018.
#'
#' @param spec a (named) list of signal-to-noise ratio data sets: each data set
#'   itself should be list containing at least a named element \code{snr} which
#'   is an object of class \code{"spec"} providing signal-to-noise ratios as a
#'   function of frequency. For Figure 3 in Münch and Laepple (2018) set
#'   \code{spec} to the output from \code{\link{PublicationSNR}}.
#' @param names an optional character vector of names of the proxy data
#'   sets. If \code{NULL}, the names of \code{spec} are used or, if not present,
#'   default names.
#' @param f.cut Shall the spectra be cut at the cutoff frequency constrained
#'   by the diffusion correction strength? Defaults to \code{FALSE}.
#' @param col a numeric or character vector of colors to use for the plotting
#'   with length recycled to match \code{length(spec)}.
#' @param plt.ann if \code{"default"} use axis annotation as in Fig. 3 of Münch
#'   and Laepple (2018). Since no other fixed annotation scheme is implemented,
#'   setting \code{plt.ann} to a different value will result in an error.
#' @inheritParams PlotArraySpectra
#'
#' @author Thomas Münch
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @examples
#'
#' # Plot Figure 3 in Münch and Laepple (2018)
#' # (DML and WAIS oxygen isotope data sets):
#'
#' # Load main spectral results
#' DWS <- WrapSpectralResults(
#'        dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
#'        diffusion = diffusion.tf,
#'        time.uncertainty = time.uncertainty.tf,
#'        df.log = c(0.15, 0.15, 0.1))
#'
#' # Calculate the final signal-to-noise ratio spectra
#' SNR <- proxysnr:::PublicationSNR(DWS)
#'
#' # Plot it
#' PlotSNR(SNR, f.cut = TRUE,
#'         names = c("DML", "WAIS"),
#'         col = c("black", "dodgerblue4"))
#'
#' @export
#'
PlotSNR <- function(spec, f.cut = FALSE,
                    names = NULL, col = 1 : length(spec),
                    xlim = c(500, 2), ylim = c(0.05, 5),
                    plt.ann = "default",
                    xlab = NULL, ylab = NULL,
                    xtm = NULL, ytm = NULL,
                    xtl = NULL, ytl = NULL) {

  if (length(col) != length(spec)) col <- rep(col, length.out = length(spec))

  # Axis settings
  
  if (plt.ann == "default") {
    set.xlab <- "Time period (yr)"
    set.ylab <- "Signal-to-Noise Ratio"
    set.xtm <- NULL
    set.xtl <- NULL
    set.ytm <- c(0.05, 0.1, 0.5, 1, 5)
    set.ytl <- set.ytm
  } else {
    stop("Invalid axis annotation setting.")
  }
  if (!is.null(xlab)) set.xlab = xlab
  if (!is.null(ylab)) set.ylab = ylab
  if (!is.null(xtm)) {set.xtm = xtm; set.xtl = xtl}
  if (!is.null(ytm)) {set.ytm = ytm; set.ytl = ytl}

  op <- graphics::par(mar = c(5, 6, 0.5, 0.5), las = 1,
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
  on.exit(graphics::par(op))

  # Plot SNR

  plot.snr <- function(snr, xlim, ylim, lwd, col,
                       conf, removeF, removeL, add = FALSE) {

    if (!add) {
      LPlot(snr, bPeriod = TRUE, bNoPlot = TRUE, axes = FALSE,
            xlim = xlim, ylim = ylim, xlab = "", ylab = "")
    }

    LLines(snr, conf = conf, bPeriod = TRUE,
           removeFirst = removeF, removeLast = removeL,
           col = col, lwd = lwd)

  }

  for (i in 1 : length(spec)) {

    add <- ifelse(i == 1, FALSE, TRUE)

    removeL <- 0
    if (f.cut) {
      idx <- spec[[i]]$f.cutoff[1]
      if (is.null(idx)) {
        warning(
          "f.cut = TRUE but no cutoff frequency specified in input.")
      } else {
        removeLast <-
          length(idx : length(spec[[i]]$snr$freq))
      }
    }
    
    plot.snr(spec[[i]]$snr, xlim = xlim, ylim = ylim, lwd = 2, col = col[i],
             conf = FALSE, removeF = 1, removeL = removeLast, add = add)

  }

  # Axis and legends settings
  
  graphics::axis(1, at = set.xtm, labels = set.xtl)
  graphics::axis(2, at = set.ytm, labels = set.ytl)
  graphics::mtext(set.xlab, side = 1, line = 3.5,
                  cex = graphics::par()$cex.lab)
  graphics::mtext(set.ylab, side = 2, line = 4.25,
                  cex = graphics::par()$cex.lab, las = 0)

  if (is.null(names)) {
    names <- names(spec)
    if (is.null(names))
      names <- paste("data", 1 : length(spec), sep = "")
  }
  if (length(names) != length(spec))
    warning("Number of data sets does not match given nuber of names.",
            call. = FALSE)
  graphics::legend("topleft", legend = names, col = col,
                   seg.len = 3, lty = 1, lwd = 2, bty = "n")
  
}
    
#' Plot proxy stack correlation
#'
#' Plot the correlation of the spatial average of a certain number of proxy
#' records with the underlying common signal depending on the number of records
#' averaged and their temporal resolution, such as in Fig. 4 in Münch and
#' Laepple, 2018.
#'
#' @param data a list of the correlation data (e.g. as output by
#'   \code{\link{ObtainStackCorrelation}}), which must have two elements:
#'   \code{freq} and \code{correlation}, where \code{freq} contains the
#'   frequency axis of the underlying proxy data set (to obtain an axis for the
#'   temporal averaging period), and where \code{correlation} is a \code{n * m}
#'   matrix of the correlation values. The number of columns of
#'   \code{correlation} must match the number of frequency values (i.e. the
#'   length of \code{freq}), and the row index stands for the number of proxy
#'   records averaged.
#' @param col.pal a color palette function to be used to assign colors in the
#'   plot; the default `NULL` means to calculate the palette function internally
#'   from ten colours of the diverging \code{RdYlBu} palette in the ColorBrewer
#'   2.0 collection.
#' @param label an optional label of the data set to be displayed at the top
#'   of the plot.
#' @param plt.ann if \code{"default"} use axis annotation as in Fig. 4 of Münch
#'   and Laepple (2018). Since no other fixed annotation scheme is implemented,
#'   setting \code{plt.ann} to a different value will result in an error.
#' @inheritParams PlotArraySpectra
#' @param xtm.min if not \code{NULL} use these specific x axis minor tick marks
#'   to override the default setting. Set to \code{NA} to omit minor ticks at
#'   all.
#' @param ytm.min as \code{xtm.min} for minor y axis tick marks.
#' @param xlim the x limits (x1, x2) of the plot. Set to \code{NA} to use
#'   default limits, or supply a numeric vector of length 2 with custom
#'   limits in log units. In the latter case, setting either of the elements to
#'   \code{NA} results in using the default limit for this element only.
#' @param ylim as \code{xlim} for the y limits of the plot.
#'
#' @author Thomas Münch
#' @seealso \code{\link{ObtainStackCorrelation}}; <https://colorbrewer2.org/>
#'
#' @references
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @examples
#'
#' # Plot Figure 5 in Münch and Laepple (2018)
#' # (DML and WAIS oxygen isotope data sets):
#'
#' # Load main spectral results
#' DWS <- WrapSpectralResults(
#'        dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
#'        diffusion = diffusion.tf,
#'        time.uncertainty = time.uncertainty.tf,
#'        df.log = c(0.15, 0.15, 0.1))
#'
#' # Calculate the final signal-to-noise ratio spectra
#' SNR <- proxysnr:::PublicationSNR(DWS)
#'
#' # Calculate the correlations
#' crl <- ObtainStackCorrelation(SNR$dml, N = 1 : 20,
#'                               freq.cut.lower = 1 / 100,
#'                               freq.cut.upper = SNR$dml$f.cutoff[2])
#'
#' # Plot it
#' PlotStackCorrelation(data = crl, label = "DML", ylim = c(NA, log(50)))
#'
#' @export
#'
PlotStackCorrelation <- function(data, col.pal = NULL,
                                 label = "",
                                 plt.ann = "default",
                                 xlab = NULL, ylab = NULL,
                                 xtm = NULL, ytm = NULL,
                                 xtl = NULL, ytl = NULL,
                                 xtm.min = NULL, ytm.min = NULL,
                                 xlim = NA, ylim = NA) {

  # Error checking
  if (!is.list(data)) stop("'data' must be a list.", call. = FALSE)
  if (!all(utils::hasName(data, c("freq", "correlation"))))
    stop("'data' must have elements 'freq' and 'correlation'.", call. = FALSE)
  if (!is.matrix(data$correlation))
    stop("Element 'correlation' must be a n x m matrix.", call. = FALSE)
  if (length(data$freq) != ncol(data$correlation)) {
    stop("Length of element 'freq' must match number of columns ",
         "of element 'correlation'.", call. = FALSE)
  }

  if (plt.ann == "default") {
    set.xlab <- "Number of cores"
    set.ylab <- "Averaging period (yr)"
    set.xtm <- c(1, 2, 5, 10, 20)
    set.xtl <- set.xtm
    set.ytm <- c(2, 5, 10, 20, 50)
    set.ytl <- set.ytm
    set.xtm.min <- c(3, 4, 6 : 9, 11 : 19)
    set.ytm.min <- c(3, 4, 6 : 9, 30, 40)
  } else {
    stop("Invalid axis annotation setting.")
  }
  if (!is.null(xlab)) set.xlab = xlab
  if (!is.null(ylab)) set.ylab = ylab
  if (!is.null(xtm)) {set.xtm = xtm; set.xtl = xtl}
  if (!is.null(ytm)) {set.ytm = ytm; set.ytl = ytl}
  if (!is.null(xtm.min)) set.xtm.min <- xtm.min
  if (!is.null(ytm.min)) set.ytm.min <- ytm.min

  # Set default color palette function

  if (!length(col.pal)) {
    col.pal <- grDevices::colorRampPalette(
                            rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
  }

  # Gather input data and transform to log scale

  n <- nrow(data$correlation)

  x <- 1 : n
  x <- log(x)
  
  y <- 2 * data$freq
  y <- rev(1 / y)
  y <- log(y)

  z <- data$correlation

  # Graphics settings
  
  op <- graphics::par(mar = c(5, 5, 2, 0.5), las = 1,
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
  on.exit(graphics::par(op))

  if (length(xlim) == 1) {
    if (is.na(xlim)) {
      xlim <- range(x, finite = TRUE)
    } else {
      stop("Invalid x limit setting.")
    }
  } else {
    idx <- which(is.na(xlim))
    xlim[idx] <- range(x, finite = TRUE)[idx]
  }
  if (length(ylim) == 1) {
    if (is.na(ylim)) {
      ylim <- range(y, finite = TRUE)
    } else {
      stop("Invalid y limit setting.")
    }
  } else {
    idx <- which(is.na(ylim))
    ylim[idx] <- range(y, finite = TRUE)[idx]
  }
  
  # Plot filled contour map
  graphics::filled.contour(x, y, z,
                           color.palette = col.pal,
                           xlim = xlim, ylim = ylim,
                           zlim = c(0, 1),
                           plot.title = graphics::title(xlab = set.xlab, ylab = set.ylab),
                           plot.axes =
                             {
                               graphics::contour(x, y, z,
                                                 add = TRUE, labcex = 1, lwd = 1);
                               graphics::axis(1, at = log(set.xtm), label = set.xtl);
                               graphics::axis(1, at = log(set.xtm.min), label = FALSE,
                                              tcl = 0.5 * graphics::par("tcl"));
                               graphics::axis(2, at = log(set.ytm), label = set.ytl);
                               graphics::axis(2, at = log(set.ytm.min), label = FALSE,
                                              tcl = 0.5 * graphics::par("tcl"));
                             }
                           )

  # Plot labels
  op.usr <- graphics::par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  graphics::text(0.98, 0.5, labels = "Correlation",
                 srt = -90, xpd = NA, cex = graphics::par()$cex.lab)
  graphics::text(0.01, 1.04, adj = c(0, 0.5), labels = label,
                 xpd = NA, cex = graphics::par()$cex.lab)
  graphics::par(op.usr)
  
}

#' Plot transfer functions
#'
#' Plot the spectral transfer functions of the effects of diffusion and time
#' uncertainty, such as in Fig. B1 in Münch and Laepple, 2018.
#'
#' @param dtf A list of transfer function data sets: each data set is an object
#'   of class \code{"spec"} (see \code{?spectrum}) with minimum components
#'   \code{freq} and \code{spec}, or simply a named list with these two, where
#'   component \code{freq} is a numeric vector providing a frequency axis and
#'   component \code{spec} a numeric vector with the corresponding diffusion
#'   transfer function values. If \code{NULL} (the default), the diffusion
#'   transfer function provided with the package is plotted, which corresponds
#'   to Figure B1 in Münch and Laepple (2018).
#' @param ttf As \code{dtf} but providing time uncertainty transfer
#'   functions. If \code{NULL} (the default), the time uncertainty transfer
#'   function provided with the package is plotted, which corresponds to Figure
#'   B1 in Münch and Laepple (2018).
#' @param names an optional character vector of names for the transfer function
#'   data sets. If \code{NULL}, the names of \code{dtf} and \code{ttf} are used
#'   or, if not present, default names. If the diffusion and time uncertainty
#'   data sets differ in number, you can provide a list of two vectors of names.
#' @param col a numeric or character vector of colors to use for the plotting;
#'   if \code{NULL} default colors are used.
#' @param dtf.threshold optional critical diffusion transfer function
#'   value to plot a corresponding horizontal line and vertical lines of
#'   corresponding frequency cutoff values (omitted for \code{NULL}).
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim1 the y limits (y1, y2) of the diffusion transfer function plot.
#' @param ylim2 the y limits (y1, y2) of the time uncertainty transfer function
#'   plot.
#' @param plt.ann if \code{"default"} use axis annotation as in Fig. B1 of
#'   Münch and Laepple (2018). Since no other fixed annotation scheme is
#'   implemented, setting \code{plt.ann} to a different value will result in an
#'   error.
#' @inheritParams PlotArraySpectra
#' @param ylab1 if not ‘NULL’ use this specific y axis label for the first plot
#'   to override the default setting.
#' @param ylab2 if not ‘NULL’ use this specific y axis label for the second
#'   plot to override the default setting.
#' @param ytm1 as \code{xtm} for the first y axis.
#' @param ytm2 as \code{xtm} for the second y axis.
#' @param ytl1 as \code{xtl} for the first y axis.
#' @param ytl2 as \code{xtl} for the second y axis.
#'
#' @author Thomas Münch
#'
#' @references Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @examples
#'
#' # Plot Figure B1 in Münch and Laepple (2018), i.e. the used transfer
#' # functions to correct the DML and WAIS isotope spectra:
#'
#' PlotTF(names = c("DML1", "DML2", "WAIS"), dtf.threshold = 0.5,
#'        col = c("black", "firebrick", "dodgerblue"))
#'
#' @export
#'
PlotTF <- function(dtf = NULL, ttf = NULL,
                   names = NULL, col = NULL,
                   dtf.threshold = NULL,
                   xlim = c(500, 2),
                   ylim1 = c(0.005, 5), ylim2 = c(0.2, 1.5),
                   plt.ann = "default",
                   xlab = NULL, ylab1 = NULL, ylab2 = NULL,
                   xtm = NULL, ytm1 = NULL, ytm2 = NULL,
                   xtl = NULL, ytl1 = NULL, ytl2 = NULL) {

  # Gather or load transfer functions
  
  if (is.null(dtf)) dtf <- proxysnr::diffusion.tf
  if (is.null(ttf)) ttf <- proxysnr::time.uncertainty.tf

  # Axis settings

  if (plt.ann == "default") {
    set.xlab <- "Time period (yr)"
    set.ylab1 <- expression(bar(G))
    set.ylab2 <- expression(Phi)
    set.xtm <- NULL
    set.xtl <- NULL
    set.ytm1 <- c(0.01, 0.05, 0.1, 0.5, 1, 5)
    set.ytl1 <- set.ytm1
    set.ytm2 <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
    set.ytl2 <- NULL
  } else {
    stop("Invalid axis annotation setting.")
  }
  if (!is.null(xlab)) set.xlab <- xlab
  if (!is.null(ylab1)) set.ylab1 <- ylab1
  if (!is.null(ylab2)) set.ylab2 <- ylab2
  if (!is.null(xtm)) {set.xtm = xtm; set.xtl = xtl}
  if (!is.null(ytm1)) {set.ytm1 = ytm1; set.ytl1 = ytl1}
  if (!is.null(ytm2)) {set.ytm2 = ytm2; set.ytl2 = ytl2}

  # Plot parameters

  op <- graphics::par(mar = c(0, 0, 0, 0), las = 1,
                      oma = c(5, 5, 0.5, 0.5), mfrow = c(2,1),
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
  on.exit(graphics::par(op))

  if (is.null(col)) col <- 1 : max(length(dtf), length(ttf))

  if (is.null(names)) {
    nam1 <- names(dtf)
    nam2 <- names(ttf)
    if (is.null(nam1))
      nam1 <- paste("data", 1 : length(dtf), sep = "")
    if (is.null(nam2))
      nam2 <- paste("data", 1 : length(dtf), sep = "")
  } else if (is.list(names)) {
    nam1 <- names[[1]]
    nam2 <- names[[2]]
  } else {
    nam1 <- nam2 <- names
  }

  if (length(nam1) != length(dtf))
    warning("dtf: Number of data sets does not match number of names.",
            call. = FALSE)
  if (length(nam2) != length(ttf))
    warning("ttf: Number of data sets does not match number of names.",
            call. = FALSE)

  if (!is.null(dtf.threshold)) {
    f.cutoff <- sapply(dtf, function(x) {
      x$freq[which(x$spec <= dtf.threshold)[1]]}
      )
  }
  
  # Wrapper function for the legend
  
  leg <- function(names, col) {
    graphics::legend("bottomleft", legend = names,
                     lwd = 2, lty = 1, col = col, bty = "n")
  }

  
  # Plot diffusion transfer functions

  ii <- length(dtf)
  jj <- length(ttf)
  
  LPlot(dtf[[1]], bNoPlot = TRUE, bPeriod = TRUE, axes = FALSE,
        xlab = "", ylab = "", xlim = xlim, ylim = ylim1)

  for (i in 1 : ii) {
    
    LLines(dtf[[i]], bPeriod = TRUE, lwd = 2, col = col[i])

    if (!is.null(dtf.threshold)) {
      graphics::lines(x = rep(1 / f.cutoff[i], 2),
                      y = c(ylim1[1] / 10, dtf.threshold),
                      lwd = 1, lty = 2, col = col[i])
    }
  }

  if (!is.null(dtf.threshold)) {
    graphics::lines(x = c(2 * xlim[1], min(1 / f.cutoff[!is.na(f.cutoff)])),
                    y = rep(dtf.threshold, 2),
                    lwd = 1, lty = 2, col = "darkgrey")
  }

  graphics::mtext("a", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)
  graphics::mtext(set.ylab1, side = 2, line = 3.5,
                  las = 0, cex = graphics::par()$cex.lab)

  graphics::axis(2, at = set.ytm1, labels = set.ytl1)
  graphics::box()

  # Place an extra legend if different number of data sets are used
  if (ii != jj) leg(nam1, col)


  # Plot time uncertainty transfer functions

  LPlot(ttf$dml1, bNoPlot = TRUE, bPeriod = TRUE, axes = FALSE,
        xlab = "", ylab = "", xlim = xlim, ylim = ylim2)

  for (i in 1 : jj) {
    LLines(ttf[[i]], bPeriod = TRUE, lwd = 2, col = col[i])
  }

  graphics::mtext("b", side = 3, adj = 0.01, padj = 0.5,
                  line = -1, font = 2, cex = graphics::par()$cex.lab)
  graphics::mtext(set.ylab2, side = 2, line = 3.5,
                  las = 0, cex = graphics::par()$cex.lab)
  graphics::mtext(set.xlab, side = 1, line = 3.5,
                  las = 0, cex = graphics::par()$cex.lab)

  graphics::axis(1, at = set.xtm, labels = set.xtl)
  graphics::axis(2, at = set.ytm2, labels = set.ytl2)
  graphics::box()

  leg(nam2, col)

}

