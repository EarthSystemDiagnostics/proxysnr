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

#' Plot spectra from spatial proxy data array
#'
#' Plot the spectral estimates from a spatial proxy data array (e.g., as in the
#' firn core analysis of Münch and Laepple, 2018, Fig. 1).
#'
#' @param spec output from \code{\link{ObtainArraySpectra}}.
#' @param marker vector of optional frequency values to mark certain parts of
#'   the plot, e.g. a cutoff frequency, in form of vertical grey dashed lines.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim the y limits (y1, y2) of the plot.
#' @param col a three-element vector with the colours for the individual
#'   spectra (\code{col[1]}), the mean (\code{col[2]}) and the stack spectrum
#'   (\code{col[3]}).
#' @param col.sn a two-element vector with the colours for the approximate
#'   signal (\code{col[1]}) and noise shading (\code{col[2]}).
#' @param alpha.sn opacity factor for the colours in \code{col.sn} within
#'   [0,1].
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param xtm x axis tick mark positions; if \code{NULL} computed by
#'   \code{\link[graphics]{axis}}.
#' @param ytm y axis tick mark positions; if \code{NULL} computed from the
#'   \code{ylim} range in steps of powers of 10.
#' @param xtl x axis tick mark labels; if \code{NULL} determined automatically
#'   from \code{xtm}, else it must be a vector of labels of the same length as
#'   \code{xtm}.
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
#' library(magrittr)
#'
#' # create artificial proxy data as a showcase dataset
#' nc <- 5
#' nt <- 1000
#' clim <- as.numeric(stats::arima.sim(model = list(ar = 0.7), n = nt))
#' noise <- as.data.frame(replicate(nc, rnorm(n = nt)))
#'
#' # array of five "cores" recording the same climate but independent noise:
#' cores <- clim + noise
#'
#' # obtain the spectra and plot them
#' cores %>%
#'   ObtainArraySpectra(df.log = 0.05) %>%
#'   PlotArraySpectra(xlim = c(500, 2))
#'
#' @export
#'
PlotArraySpectra <- function(spec, marker = NA,
                             xlim = c(100, 2), ylim = c(0.005, 50),
                             col = c("darkgrey", "black", "burlywood4"),
                             col.sn = c("dodgerblue4", "firebrick4"),
                             alpha.sn = 0.2,
                             xlab = "Time period (yr)",
                             ylab = "Power spectral density",
                             xtm = NULL, ytm = NULL,
                             xtl = NULL, ytl = NULL) {

  # Error checking

  if (!is.list(spec)) stop("`spec` must be a list.", call. = FALSE)
  if (!all(utils::hasName(spec, c("single", "mean", "stack")))) {
      stop("`spec` must have elements `single`, `mean` and `stack`.",
           call. = FALSE)
  }
  is.spectrum(spec$mean)
  is.spectrum(spec$stack)
  if (!is.list(spec$single))
    stop("`spec$single` must be a list of spectra.", call. = FALSE)

  # Gather input
  
  N <- length(spec$single)    
  psd1 <- spec$mean
  psd2 <- spec$stack
  psd3 <- psd1
  psd3$spec <- psd3$spec / N

  # Axis settings

  if (is.null(xtl)) xtl <- xtm
  if (is.null(ytm)) ytm <- 10^(seq(log10(ylim[1]), log10(ylim[2]), by = 1))
  if (is.null(ytl))
    ytl <- format(ytm, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  
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

  # Horizontal marker line(s)

  lapply(marker, function(m) {
    graphics::lines(x = rep(1 / m, 2), y = c(ylim[1]/10, ylim[2]),
                    lty = 2, col = "darkgrey")
  })

  # Individual spectra

  for (i in 1 : N) {

    tryCatch(
    {
      LLines(spec$single[[i]], conf = FALSE, bPeriod = TRUE,
             removeFirst = 1, removeLast = 1,
             col = col[1])
    }, error = function(cond) {
      msg <- sprintf("Cannot plot `spec$single[[%s]]`: no spectral object.", i)
      stop(msg, call. = FALSE)
    })
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

  graphics::axis(1, at = xtm, labels = xtl)
  graphics::axis(2, at = ytm, labels = ytl)

  graphics::mtext(xlab, side = 1, line = 3.5, cex = graphics::par()$cex.lab)
  graphics::mtext(ylab, side = 2, line = 4.5, cex = graphics::par()$cex.lab,
                  las = 0)

  graphics::legend("bottomleft",
                   c("Individual spectra", "Mean spectrum",
                     "Spectrum of stacked record", "Mean spectrum scaled by 1/n"),
                   col = c(col[1 : 3], col[2]),
                   lty = c(1, 1, 1, 5),
                   lwd = c(1, 2, 2, 1), seg.len = 2.5, bty = "n")

}

#' Plot proxy signal-to-noise ratios
#'
#' Plot proxy signal-to-noise ratios of several datasets as a function of
#' timescale (e.g., as in the firn core analysis of Münch and Laepple, 2018,
#' Fig. 3).
#'
#' @param spec a (named) list of signal-to-noise ratio data sets: each data set
#'   itself should be list containing at least a named element \code{snr} which
#'   is an object of class \code{"spec"} providing signal-to-noise ratios as a
#'   function of frequency. For Figure 3 in Münch and Laepple (2018) set
#'   \code{spec} to the output from \code{\link{PublicationSNR}}.
#' @param f.cut Shall the spectra be cut at the cutoff frequency constrained
#'   by the diffusion correction strength? Defaults to \code{FALSE}.
#' @param names an optional character vector of names of the proxy data
#'   sets. If \code{NULL}, the names of \code{spec} are used or, if not present,
#'   default names.
#' @param col a numeric or character vector of colors to use for the plotting
#'   with length recycled to match \code{length(spec)}.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param ytm y axis tick mark positions; default setting (\code{NULL}) uses
#'   \code{c(0.05, 0.1, 0.5, 1, 5)}.
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
#' # create toy data
#' n <- 100
#' spec <- list(
#'   data1 = list(snr = list(freq = seq(0.01, 0.5, length.out = n),
#'                           spec = seq(1, 0.1, length.out = n))),
#'   data2 = list(snr = list(freq = seq(0.005, 0.5, length.out = n),
#'                           spec = seq(5, 0.1, length.out = n)))
#' )
#'
#' # plot SNR data
#' PlotSNR(spec)
#'
#' @export
#'
PlotSNR <- function(spec, f.cut = FALSE,
                    names = NULL, col = 1 : length(spec),
                    xlim = c(500, 2), ylim = c(0.05, 5),
                    xlab = "Time period (yr)", ylab = "Signal-to-Noise Ratio",
                    xtm = NULL, ytm = NULL,
                    xtl = NULL, ytl = NULL) {

  # Error checking
  if (!is.list(spec)) stop("`spec` must be a list.", call. = FALSE)

  if (length(col) != length(spec)) col <- rep(col, length.out = length(spec))

  # Graphics settings

  if (is.null(xtl)) xtl <- xtm
  if (is.null(ytm)) ytm <- c(0.05, 0.1, 0.5, 1, 5)
  if (is.null(ytl)) ytl <- ytm

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

  if (f.cut) {

    if (any(sapply(spec, function(x) {!utils::hasName(x, "f.cutoff")})))
      warning("`f.cut = TRUE` but cutoff frequency missing in input.",
              call. = FALSE)
  }

  for (i in 1 : length(spec)) {

    add <- ifelse(i == 1, FALSE, TRUE)

    removeLast <- 0
    if (f.cut) {
      idx <- spec[[i]]$f.cutoff[1]
      if (!is.null(idx)) {
        removeLast <- length(idx : length(spec[[i]]$snr$freq))
      }
    }

    if (!is.list(spec[[i]]))
      stop(sprintf("`spec[[%s]]` must be a list.", i), call. = FALSE)

    if (!utils::hasName(spec[[i]], "snr")) {
      stop(sprintf("`spec[[%s]]` must have an element `snr`.", i),
           call. = FALSE)
    }

    tryCatch(is.spectrum(spec[[i]]$snr), error = function(cond) {
      stop(sprintf("Cannot plot `spec[[%s]]$snr`: no spectral object.", i),
           call. = FALSE)
    })
    
    plot.snr(spec[[i]]$snr, xlim = xlim, ylim = ylim, lwd = 2, col = col[i],
             conf = FALSE, removeF = 1, removeL = removeLast, add = add)

  }

  # Axis and legends settings
  
  graphics::axis(1, at = xtm, labels = xtl)
  graphics::axis(2, at = ytm, labels = ytl)
  graphics::mtext(xlab, side = 1, line = 3.5, cex = graphics::par()$cex.lab)
  graphics::mtext(ylab, side = 2, line = 4.25, cex = graphics::par()$cex.lab,
                  las = 0)

  if (is.null(names)) {
    names <- names(spec)
    if (is.null(names))
      names <- paste("data", 1 : length(spec), sep = "")
  }
  if (length(names) != length(spec))
    warning("Number of data sets does not match supplied number of names.",
            call. = FALSE)
  graphics::legend("topleft", legend = names, col = col,
                   seg.len = 3, lty = 1, lwd = 2, bty = "n")
  
}
    
#' Plot proxy stack correlation
#'
#' Plot the correlation of the spatial average of a certain number of proxy
#' records with the underlying common signal depending on the number of records
#' averaged and their temporal resolution (e.g., as in the firn core analysis of
#' Münch and Laepple, 2018, Fig. 4).
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
#' @param xlim the x limits (x1, x2) of the plot. Set to \code{NA} to use
#'   default limits calculated from the x data range, or supply a numeric vector
#'   of length 2 with custom limits. In the latter case, setting either of the
#'   elements to \code{NA} results in using the default limit for this element
#'   only; see the example.
#' @param ylim as \code{xlim} for the y limits of the plot.
#' @param xtm x axis tick mark positions; default setting (\code{NULL}) uses
#'   \code{c(1, 2, 5, 10, 20)}.
#' @param ytm y axis tick mark positions; default setting (\code{NULL}) uses
#'   \code{c(2, 5, 10, 20, 50)}.
#' @inheritParams PlotSNR
#' @param xtm.min x axis minor tick marks; default setting \code{NULL} uses
#'   minor tick marks that are adapted to the default x major tick marks. Set to
#'   \code{NA} to omit minor ticks at all.
#' @param ytm.min as \code{xtm.min} for minor y axis tick marks.
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
#' # create a toy correlation dataset, which mimicks an increase
#' # in correlation with timescale and with the number of cores
#' # averaged, and plot it:
#'
#' nf <- 20
#' nc <- 5
#'
#' data <- list(
#'   freq = seq(0.01, 0.5, length.out = nf),
#'   correlation =
#'     matrix(seq(0.05, 0.9, length.out = nf), nrow = nc, ncol = nf,
#'            byrow = TRUE) +
#'     matrix(seq(0, 0.3, length.out = nc), nrow = nc, ncol = nf) +
#'     matrix(rnorm(nf * nc, sd = 0.02), nrow = nc, ncol = nf)
#' )
#' data$correlation[data$correlation > 1] <- 1
#' data$correlation[data$correlation < 0] <- 0
#'
#' PlotStackCorrelation(data)
#'
#'
#' @export
#'
PlotStackCorrelation <- function(data, col.pal = NULL, label = "",
                                 xlim = NA, ylim = NA,
                                 xlab = "Number of cores",
                                 ylab = "Averaging period (yr)",
                                 xtm = NULL, ytm = NULL,
                                 xtl = NULL, ytl = NULL,
                                 xtm.min = NULL, ytm.min = NULL) {

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

  if (is.null(xtm)) xtm <- c(1, 2, 5, 10, 20)
  if (is.null(ytm)) ytm <- c(2, 5, 10, 20, 50)
  if (is.null(xtl)) xtl <- xtm
  if (is.null(ytl)) ytl <- ytm
  if (is.null(xtm.min)) xtm.min <- c(3, 4, 6 : 9, 11 : 19)
  if (is.null(ytm.min)) ytm.min <- c(3, 4, 6 : 9, 30, 40)

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
    xlim <- log(xlim)
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
    ylim <- log(ylim)
    idx <- which(is.na(ylim))
    ylim[idx] <- range(y, finite = TRUE)[idx]
  }
  
  # Plot filled contour map
  graphics::filled.contour(x, y, z,
                           color.palette = col.pal,
                           xlim = xlim, ylim = ylim,
                           zlim = c(0, 1),
                           plot.title = graphics::title(xlab = xlab, ylab = ylab),
                           plot.axes =
                             {
                               graphics::contour(x, y, z,
                                                 add = TRUE, labcex = 1, lwd = 1);
                               graphics::axis(1, at = log(xtm), label = xtl);
                               graphics::axis(1, at = log(xtm.min), label = FALSE,
                                              tcl = 0.5 * graphics::par("tcl"));
                               graphics::axis(2, at = log(ytm), label = ytl);
                               graphics::axis(2, at = log(ytm.min), label = FALSE,
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
#' uncertainty (e.g., as in the analysis of Münch and Laepple, 2018, Fig. B1).
#'
#' @param dtf A list of transfer function data sets: each data set is an object
#'   of class \code{"spec"} (see \code{?spectrum}) with minimum components
#'   \code{freq} and \code{spec}, or simply a named list with these two, where
#'   component \code{freq} is a numeric vector providing a frequency axis and
#'   component \code{spec} a numeric vector with the corresponding diffusion
#'   transfer function values. If \code{NULL} (the default), the diffusion
#'   transfer functions provided with the package are plotted, which corresponds
#'   to Figure B1 in Münch and Laepple (2018).
#' @param ttf As \code{dtf} but providing time uncertainty transfer
#'   functions. If \code{NULL} (the default), the time uncertainty transfer
#'   functions provided with the package are plotted, which corresponds to
#'   Figure B1 in Münch and Laepple (2018).
#' @param names an optional character vector of names for the transfer function
#'   data sets. If \code{NULL}, the names of \code{dtf} and \code{ttf} are used
#'   or, if not present, default names. If the diffusion and time uncertainty
#'   data sets differ in number, provide a list of two vectors of names.
#' @param col a numeric or character vector of colors to use for the plotting;
#'   if \code{NULL} default colors are used.
#' @param dtf.threshold optional critical diffusion transfer function
#'   value to plot a corresponding horizontal line and vertical lines of
#'   corresponding frequency cutoff values (omitted for \code{NULL}).
#' @param xlab x axis label.
#' @param ylab1 y axis label for the diffusion transfer function.
#' @param ylab2 y axis label for the time uncertainty transfer function.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim1 the y limits (y1, y2) of the diffusion transfer function plot.
#' @param ylim2 the y limits (y1, y2) of the time uncertainty transfer function
#'   plot.
#' @param ytm1 y axis tick mark positions on the diffusion transfer function
#'   plot; default setting (\code{NULL}) uses
#'   \code{c(0.01, 0.05, 0.1, 0.5, 1, 5)}.
#' @param ytm2 y axis tick mark positions on the time uncertainty transfer
#'   function plot; default setting (\code{NULL}) uses
#'   \code{c(0.2, 0.4, 0.6, 0.8, 1, 1.2)}.
#' @inheritParams PlotArraySpectra
#' @param ytl1 y axis tick mark labels on the diffusion transfer function plot;
#'   if \code{NULL} determined automatically from \code{ytm1}, else it must be a
#'   vector of labels of the same length as \code{ytm1}.
#' @param ytl2 as \code{ytl1} for the y axis on the time uncertainty transfer
#'   function plot.
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
                   xlab = "Time period (yr)",
                   ylab1 = expression(bar(G)), ylab2 = expression(Phi),
                   xlim = c(500, 2), ylim1 = c(0.005, 5), ylim2 = c(0.2, 1.5),
                   xtm = NULL, ytm1 = NULL, ytm2 = NULL,
                   xtl = NULL, ytl1 = NULL, ytl2 = NULL) {

  # Gather or load transfer functions

  if (is.null(dtf) & is.null(ttf)) {
    dtf <- proxysnr::diffusion.tf
    ttf <- proxysnr::time.uncertainty.tf
  }

  plot.dtf <- !is.null(dtf)
  plot.ttf <- !is.null(ttf)

  plot.both <- plot.dtf & plot.ttf

  # Axis settings

  if (is.null(xtl)) xtl <- xtm
  if (is.null(ytm1)) ytm1 <- c(0.01, 0.05, 0.1, 0.5, 1, 5)
  if (is.null(ytl1)) ytl1 <- ytm1
  if (is.null(ytm2)) ytm2 <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
  if (is.null(ytl1)) ytl1 <- ytm1

  # Plot parameters

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

  # wrapper function for plot legend

  leg <- function(names, col) {
    graphics::legend("bottomleft", legend = names,
                     lwd = 2, lty = 1, col = col, bty = "n")
  }

  # wrapper function to plot transfer functions

  .plottf <- function(tf, col, xlab, ylab, xlim, ylim, xtm, ytm, xtl, ytl,
                      nam, dtf.threshold = NULL, plot.legend = TRUE,
                      plot.xax = TRUE, ch = "") {

    if (length(nam) != length(tf)) {
      warning(
        sprintf("%s: Number of data sets does not match number of names.", ch),
        call. = FALSE)
    }

    LPlot(tf[[1]], bNoPlot = TRUE, bPeriod = TRUE, axes = FALSE,
          xlab = "", ylab = "", xlim = xlim, ylim = ylim)

    sapply(seq(length(tf)), function(i) {
      LLines(tf[[i]], bPeriod = TRUE, lwd = 2, col = col[i])
    })

    if (!is.null(dtf.threshold)) {

      f.cutoff <- sapply(tf, function(x) {
        x$freq[which(x$spec <= dtf.threshold)[1]]
      })

      sapply(seq(length(tf)), function(i) {
        graphics::lines(x = rep(1 / f.cutoff[i], 2),
                        y = c(ylim[1] / 10, dtf.threshold),
                        lwd = 1, lty = 2, col = col[i])
      })

      graphics::lines(x = c(2 * xlim[1], min(1 / f.cutoff[!is.na(f.cutoff)])),
                      y = rep(dtf.threshold, 2),
                      lwd = 1, lty = 2, col = "darkgrey")

    }

    if (plot.xax) {

      graphics::axis(1, at = xtm, labels = xtl)
      graphics::mtext(xlab, side = 1, line = 3.5, las = 0,
                      cex = graphics::par()$cex.lab)

    }

    graphics::axis(2, at = ytm, labels = ytl)
    graphics::mtext(ylab, side = 2, line = 3.5, las = 0,
                    cex = graphics::par()$cex.lab)

    graphics::box()

    if (plot.legend) leg(nam, col)

  }

  if (!plot.both) {

    op <- graphics::par(mar = c(0, 0, 0, 0), las = 1,
                        oma = c(5, 5, 0.5, 0.5),
                        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
    on.exit(graphics::par(op))

    if (plot.dtf) {

      # Plot diffusion transfer functions

      .plottf(dtf, col = col, xlab = xlab, ylab = ylab1,
              xlim = xlim, ylim = ylim1, xtm = xtm, ytm = ytm1,
              xtl = xtl, ytl = ytl1, nam = nam1,
              dtf.threshold = dtf.threshold, plot.legend = TRUE, ch = "dtf")

    }

    if (plot.ttf) {

      # Plot time uncertainty transfer functions

      .plottf(ttf, col = col, xlab = xlab, ylab = ylab2,
              xlim = xlim, ylim = ylim2, xtm = xtm, ytm = ytm2,
              xtl = xtl, ytl = ytl2, nam = nam2,
              plot.legend = TRUE, ch = "ttf")

    }

  } else {

    # Plot both transfer functions

    op <- graphics::par(mar = c(0, 0, 0, 0), las = 1,
                        oma = c(5, 5, 0.5, 0.5), mfrow = c(2,1),
                        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25)
    on.exit(graphics::par(op))

    plot.legend <- length(dtf) != length(ttf)

    .plottf(dtf, col = col, xlab = xlab, ylab = ylab1,
            xlim = xlim, ylim = ylim1, xtm = xtm, ytm = ytm1,
            xtl = xtl, ytl = ytl1, nam = nam1,
            dtf.threshold = dtf.threshold, plot.legend = plot.legend,
            plot.xax = FALSE, ch = "dtf")

    graphics::mtext("a", side = 3, adj = 0.01, padj = 0.5,
                    line = -1, font = 2, cex = graphics::par()$cex.lab)

    .plottf(ttf, col = col, xlab = xlab, ylab = ylab2,
            xlim = xlim, ylim = ylim2, xtm = xtm, ytm = ytm2,
            xtl = xtl, ytl = ytl2, nam = nam2,
            plot.legend = TRUE, ch = "ttf")

    graphics::mtext("b", side = 3, adj = 0.01, padj = 0.5,
                    line = -1, font = 2, cex = graphics::par()$cex.lab)

  }

}
