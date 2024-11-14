<!-- badges: start -->
[![Github_Status_Badge](https://img.shields.io/badge/Github-1.0.0-blue.svg)](https://github.com/EarthSystemDiagnostics/proxysnr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2027638.svg)](https://doi.org/10.5281/zenodo.2027638)
[![R-CMD-check](https://github.com/EarthSystemDiagnostics/proxysnr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EarthSystemDiagnostics/proxysnr/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/EarthSystemDiagnostics/proxysnr/branch/dev-unit-tests/graph/badge.svg)](https://codecov.io/gh/EarthSystemDiagnostics/proxysnr)
<!-- badges: end -->

# proxysnr: Separate Signal and Noise in Climate Proxy Records

## Introduction

`proxysnr` implements a method working in the spectral domain to separate the
common signal from the local noise as recorded by a spatial network of climate
proxy records, which yields an estimate of the timescale dependence of the proxy
signal-to-noise ratio (SNR). The method allows the correction of the estimated
spectra for the effects of time uncertainty, proxy smoothing processes
(e.g. diffusion), and measurement noise.

The method is introduced and in detail explained in Münch and Laepple (2018),
and it has been applied there to oxygen isotope records from Antarctic firn and
ice cores.

`proxysnr` has been implemented by Dr. Thomas Münch with contributions by
Dr. Thomas Laepple, Dr. Andrew Dolman, Dr. Torben Kunz, and Dr. Mara
McPartland. Please contact Dr. Thomas Münch <<thomas.muench@awi.de>> at the
Alfred Wegener Institute, Helmholtz Centre for Polar and Marine Research,
Germany, for further information.

This work was supported by the European Research Council (ERC) under the European
Union’s Horizon 2020 research and innovation programme (grant agreement
no. 716092) and Helmholtz funding through the Polar Regions and
Coasts in the Changing Earth System (PACES) programme of the Alfred Wegener
Institute. It further contributes to the German BMBF project PalMod.

## Installation

The current version of`proxysnr` can be installed directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("EarthSystemDiagnostics/proxysnr")
```

Released versions of `proxysnr` can be downloaded in source form from its
Zenodo repository under the DOI
[10.5281/zenodo.2027638](https://doi.org/10.5281/zenodo.2027638). The package
can then be installed using

```r
# install.packages("devtools")
devtools::install(.)
```

where `.` is your source package root directory.

## Documentation and vignettes

The full package documentation is available from the [proxysnr
homepage](https://earthsystemdiagnostics.github.io/proxysnr/).

Among this documentation, `proxysnr` includes three vignettes to highlight the
main aspects of the package:

* The `vignette("proxysnr")` introduces the main signal and noise
  decomposition method with a simple surrogate data example.
* The `vignette("plot-muench-laepple-figures")` applies this method on two
  real-world proxy datasets from Antarctic firn and ice cores, thereby showing
  how to reproduce the results of Münch and Laepple (2018).
* The `vignette("calculate-transfer-functions")` demonstrates how to obtain
   spectral transfer functions describing the loss in spectral power for two
   special cases: time-uncertainty in layer-counted chronologies and isotope
   diffusion in polar firn.

## Literature cited

Münch, T. and Laepple, T.: What climate signal is contained in decadal- to
centennial-scale isotope variations from Antarctic ice cores?, Clim. Past, 14,
2053-2070, doi:
[10.5194/cp-14-2053-2018](https://doi.org/10.5194/cp-14-2053-2018), 2018.

