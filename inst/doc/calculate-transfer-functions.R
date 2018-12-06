## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE, warning = FALSE-----------------------------------
library(proxysnr)

## ------------------------------------------------------------------------
# check if package is available
pkg.flag <- requireNamespace("simproxyage", quietly = TRUE)

if (!pkg.flag) {
   message(paste("Package \"simproxyage\" not found.",
                 "No time uncertainty transfer functions will be calculated."))
} else {
  library(simproxyage)
}


## ------------------------------------------------------------------------
set.seed(4122018)

## ------------------------------------------------------------------------
# for diffusion transfer function
dtf <- list()
# for time uncertainty transfer function
ttf <- list()

## ------------------------------------------------------------------------
res <- 1/2

## ------------------------------------------------------------------------
ff <- list()
ff$dml1 <- diffusion.tf$dml1$freq
ff$dml2 <- diffusion.tf$dml2$freq
ff$wais <- diffusion.tf$wais$freq

## ------------------------------------------------------------------------
ns <- 10

## ------------------------------------------------------------------------
nc <- 15

## ------------------------------------------------------------------------
# set time axis
t <- seq(1994, 1801)
# number of proxy data points
nt <- length(t)

## ------------------------------------------------------------------------
acp <- c(1994, 1884, 1816, 1810)

## ------------------------------------------------------------------------
rate <- 0.013

## ------------------------------------------------------------------------
tmp <- DiffusionTF(nt = nt / res, nc = nc, ns = ns,
                   sigma = diffusion.length$dml1[, -1],
                   res = res)$ratio

## ------------------------------------------------------------------------
dtf$dml1 <- tmp
dtf$dml1$freq <- ff$dml1
dtf$dml1$spec <- approx(tmp$freq, tmp$spec, dtf$dml1$freq)$y

attr(dtf$dml1, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(dtf$dml1, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(dtf$dml1, "log-smooth") <- "Log-smooth applied: No."

class(dtf$dml1) <- "spec"

## ------------------------------------------------------------------------
if (pkg.flag) {
   ttf$dml1 <- TimeUncertaintyTF(t = t, acp = acp, nt = nt, nc = nc, ns = ns,
                                 rate = rate, surrogate.fun = rnorm)$ratio
} else {
  ttf$dml1$freq <- dtf$dml1$freq
  ttf$dml1$spec <- rep(NA, length(dtf$dml1$freq))
}

attr(ttf$dml1, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(ttf$dml1, "rate") <- sprintf("Process rate used: %1.3f.", rate)
attr(ttf$dml1, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(ttf$dml1, "log-smooth") <- "Log-smooth applied: No."

class(ttf$dml1) <- "spec"

## ------------------------------------------------------------------------
nc <- 3
t <- seq(1994, 1000)
nt <- length(t)
acp <- c(1994, 1884, 1816, 1810, 1459, 1259)
rate <- 0.013

tmp <- DiffusionTF(nt = nt / res, nc = nc, ns = ns,
                   sigma = diffusion.length$dml2[, -1],
                   res = res)$ratio

dtf$dml2 <- tmp
dtf$dml2$freq <- ff$dml2
dtf$dml2$spec <- approx(tmp$freq, tmp$spec, dtf$dml2$freq)$y

if (pkg.flag) {
   ttf$dml2 <- TimeUncertaintyTF(t = t, acp = acp, nt = nt, nc = nc, ns = ns,
                                 rate = rate, surrogate.fun = rnorm)$ratio
} else {
  ttf$dml2$freq <- dtf$dml2$freq
  ttf$dml2$spec <- rep(NA, length(dtf$dml2$freq))
}

attr(dtf$dml2, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(dtf$dml2, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(dtf$dml2, "log-smooth") <- "Log-smooth applied: No."

attr(ttf$dml2, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(ttf$dml2, "rate") <- sprintf("Process rate used: %1.3f.", rate)
attr(ttf$dml2, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(ttf$dml2, "log-smooth") <- "Log-smooth applied: No."

class(dtf$dml2) <- "spec"
class(ttf$dml2) <- "spec"

## ------------------------------------------------------------------------
nc <- 5
t <- seq(2000, 1800)
nt <- length(t)
acp <- c(2000, 1992, 1964, 1885, 1837, 1816, 1810)
rate <- 0.027

tmp <- DiffusionTF(nt = nt / res, nc = nc, ns = ns,
                   sigma = diffusion.length$wais[, -1],
                   res = res)$ratio

dtf$wais <- tmp
dtf$wais$freq <- ff$wais
dtf$wais$spec <- approx(tmp$freq, tmp$spec, dtf$wais$freq)$y

if (pkg.flag) {
   ttf$wais <- TimeUncertaintyTF(t = t, acp = acp, nt = nt, nc = nc, ns = ns,
                                 rate = rate, surrogate.fun = rnorm)$ratio
} else {
  ttf$wais$freq <- dtf$wais$freq
  ttf$wais$spec <- rep(NA, length(dtf$wais$freq))
}

attr(dtf$wais, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(dtf$wais, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(dtf$wais, "log-smooth") <- "Log-smooth applied: No."

attr(ttf$wais, "version") <- sprintf("Creation date: %s.", Sys.time())
attr(ttf$wais, "rate") <- sprintf("Process rate used: %1.3f.", rate)
attr(ttf$wais, "N.sim") <- sprintf("Number of simulations used: N = %s",
	       		   	   formatC(ns, big.mark = ",", format = "d"))
attr(ttf$wais, "log-smooth") <- "Log-smooth applied: No."

class(dtf$wais) <- "spec"
class(ttf$wais) <- "spec"

## ---- fig.width = 7, fig.height = 8, fig.cap = "Estimates of the spectral transfer functions for the effects of site-specific diffusion (**a**) and time uncertainty (**b**) for three arrays of firn-core isotope records: DML1 (black), DML2 (red) and WAIS (blue) (as studied in Münch and Laepple, 2018). Plotted is in each case the average transfer function from 10 Monte Carlo simulations."----
PlotTF(dtf = dtf, ttf = ttf,
       names = c("DML1", "DML2", "WAIS"),
       col = c("black", "firebrick", "dodgerblue"))

## ---- eval = FALSE-------------------------------------------------------
#  PlotTF(names = c("DML1", "DML2", "WAIS"),
#         col = c("black", "firebrick", "dodgerblue"))

