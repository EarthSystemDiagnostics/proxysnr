# proxysnr: Separate Signal and Noise in Climate Proxy Records

------------------------------

## Introduction

**proxysnr** implements a method working in the spectral domain to separate the
common signal from the local noise as recorded by a spatial network of climate
proxy records, which yields an estimate of the timescale dependence of the proxy
signal-to-noise ratio (SNR). The implemented method includes the correction of
relevant estimated power spectral densities for the loss in spectral power by
two effects: (1) time uncertainty in the case of layer-counted record
chronologies and (2) water vapour diffusion through the open-porous firn in the
case of stable isotope records from firn and ice cores.

The method is in detail explained in Münch and Laepple (2018) and has been
applied there to Antarctic firn-core oxygen isotope records.

The R code has been implemented by Dr. Thomas Münch with contributions by
Dr. Thomas Laepple and Dr. Andrew Dolman. Please contact Dr. Thomas Münch
<<thomas.muench@awi.de>> at the Alfred Wegener Institute, Helmholtz Centre for
Polar and Marine Research, Germany, for further information.

This work was supported by Helmholtz funding through the Polar Regions and
Coasts in the Changing Earth System (PACES) programme of the Alfred Wegener
Institute, by the Initiative and Networking Fund of the Helmholtz Association
Grant VG-NH900 and by the European Research Council (ERC) under the European
Union’s Horizon 2020 research and innovation programme (grant agreement
no. 716092). It further contributes to the German BMBF project PalMod.

## Installation

**proxysnr** can be installed directly from GitHub:

```r
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("EarthSystemDiagnostics/proxysnr")
```

## Examples

Two example vignettes are provided along with the package source as rendered
`.html` files to demonstrate the main functions of the package. These vignettes
can be found in the directory `./doc/`, where `.` stands for the root of the
package source directory, while their creating R markdown source codes are
located under `./vignettes/`.

* The vignette
  [plot-muench-laepple-figures](http://htmlpreview.github.io/?https://github.com/EarthSystemDiagnostics/proxysnr/blob/master/doc/plot-muench-laepple-figures.html)
  shows the basic way of applying the package to obtain estimates of the signal,
  noise and SNR spectra. For this, the oxygen isotope data from Münch and
  Laepple (2018) are used, which are provided with this package. The vignette
  further demonstrates the plotting of all main figures shown in Münch and
  Laepple (2018).

* The vignette
   [calculate-transfer-functions](http://htmlpreview.github.io/?https://github.com/EarthSystemDiagnostics/proxysnr/blob/master/doc/calculate-transfer-functions.html)
   demonstrates the application of package functions for obtaining the spectral
   transfer functions describing the loss in spectral power by the effects of
   time uncertainty and diffusion. These transfer functions can be used to
   correct the estimated signal, noise and SNR spectra as explained in Münch and
   Laepple (2018).

After installing the package with `remotes::install_github` or after cloning
the git repository and installing the package using `devtools::install`, setting
`build_vignettes = TRUE` in both cases, the vignettes belonging to the version
of the installed package are also available directly from the `R` command line
by typing
```r
vignette("plot-muench-laepple-figures", package = "proxysnr")
```
and
```r
vignette("calculate-transfer-functions", package = "proxysnr")
```

Please note that for rebuilding the vignettes, e.g., in order to test the package
functionality, you will need to install the R package **simproxyage** used for
modelling the time uncertainty of the layer-counted isotope records in order to
calculate respective time uncertainty transfer functions. **simproxyage** is
available from [GitHub](https://github.com/EarthSystemDiagnostics/simproxyage)
and can be installed directly using
`remotes::install_github("EarthSystemDiagnostics/simproxyage")`. Also note that
the package vignettes are based on R Markdown v2 and require
[Pandoc](http://pandoc.org).

## Literature cited

Münch, T. and Laepple, T.: What climate signal is contained in decadal- to
centennial-scale isotope variations from Antarctic ice cores?, Clim. Past, 14,
2053-2070, doi:
[10.5194/cp-14-2053-2018](https://doi.org/10.5194/cp-14-2053-2018), 2018.

