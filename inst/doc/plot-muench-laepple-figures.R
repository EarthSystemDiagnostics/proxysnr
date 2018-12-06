## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE, warning = FALSE-----------------------------------
library(proxysnr)
library(RColorBrewer)

## ------------------------------------------------------------------------
DWS <- WrapSpectralResults(
    dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
    diffusion = diffusion.tf,
    time.uncertainty = time.uncertainty.tf,
    df.log = c(0.15, 0.15, 0.1))

## ------------------------------------------------------------------------
ls.str(DWS)

## ---- fig.width = 7, fig.height = 6--------------------------------------
PlotArraySpectra(ArraySpectra(dml$dml1, df.log = 0.12),
                 f.cutoff = DWS$dml1$f.cutoff[2])

## ---- warning = FALSE, message = FALSE, fig.width = 14, fig.height = 10.36----
proxysnr:::muench_laepple_fig02(DWS, f.cut = TRUE)

## ------------------------------------------------------------------------
SNR <- proxysnr:::PublicationSNR(DWS)

## ---- fig.width = 7, fig.height = 6--------------------------------------
PlotSNR(SNR, f.cut = TRUE,
        names = c("DML", "WAIS"), col = c("black", "dodgerblue4"))

## ------------------------------------------------------------------------
# for the DMl data
crl1 <- StackCorrelation(SNR$dml, N = 20,
                         freq.cut.lower = 1 / 100,
                         freq.cut.upper = SNR$dml$f.cutoff[2])

# for the WAIS data
crl2 <- StackCorrelation(SNR$wais, N = 20,
                         freq.cut.lower = 1 / 100,
                         freq.cut.upper = SNR$wais$f.cutoff[2])

## ------------------------------------------------------------------------
palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))

## ---- fig.width = 8, fig.height = 6.3------------------------------------
PlotStackCorrelation(freq = crl1$signal$freq, correlation = crl1$corr,
                     col.pal = palette, label = expression(bold("a.")~"DML"),
                     ylim = c(NA, log(50)))

## ---- fig.width = 8, fig.height = 6.3------------------------------------
PlotStackCorrelation(freq = crl2$signal$freq, correlation = crl2$corr,
                     col.pal = palette, label = expression(bold("b.")~"WAIS"),
                     ylim = c(NA, log(50)))

## ------------------------------------------------------------------------
TNS <- proxysnr:::TrenchNoise()

## ---- fig.width = 7, fig.height = 6--------------------------------------
proxysnr:::muench_laepple_fig05(SNR, TNS, f.cut = TRUE)

