# proxysnr 1.0.0

## Breaking changes

* Function names now all follow the recommended style of using verb phrases;
  hence, the following functions have been **deprecated** and replaced:

  * `ArraySpectra()` -> `ObtainArraySpectra()`
  * `SeparateSpectra()` -> `SeparateSignalFromNoise()`
  * `StackCorrelation()` -> `ObtainStackCorrelation()`
  * `IntegratedSNR()` -> `GetIntegratedSNR()`
  * `DiffusionTF()` -> `CalculateDiffusionTF()`
  * `TimeUncertaintyTF()` -> `CalculateTimeUncertaintyTF()`

  Calls that use the old name will still work, as they are piped to the respective
  function of the new name but with a deprecation warning being issued.

* Handling of transfer function data
  ([#12](https://github.com/EarthSystemDiagnostics/proxysnr/pull/12)):

  * `SeparateSignalFromNoise()`, `WrapSpectralResults()`: the input transfer
	functions now have to be supplied as spectral objects (i.e. lists with	a
	frequency and a spectral value vector of the same length);
  * following this, `WrapSpectralResults()` no longer offers the possibility of
    directly supplying the respective correction vector via signalling
    `inverse.tf = TRUE` and this parameter has been removed accordingly, and
    `SeparateSignalFromNoise()` expects a list input (i.e. the spectral object
    of the transfer function) and no longer a vector of correction (= inverse
    transfer function spectral) values;
 * to allow for more flexibility in transfer function input, the transfer
   function spectral objects are internally interpolated to the target frequency
   axis of the proxy spectral data, so e.g. the diffusion transfer function can
   be calculated on a higher temporal resolution and still directly supplied to
   `SeparateSignalFromNoise()` without requiring manual interpolation as
   previously and thus reducing code overhead for the user (see
   `vignette("calculate-transfer-functions")`. This internal interpolation of
   course requires sufficient frequency overlap between transfer function and
   proxy data, which is checked and an error being issued if it is not the case;
 * `CalculateDiffusionTF()`, `CalculateTimeUncertaintyTF()`: per default, both
   functions no longer return the entire underlying spectral data but only the
   transfer function spectral object, with the previous default return behaviour
   available via actively setting the new function parameter `verbose.output` to
   `TRUE`;
 * related function documentation was updated accordingly.

* Reorganisation of the package datasets:

  The DML and WAIS isotope datasets and the corresponding diffusion and time
  uncertainty transfer function datasets are now the only external
  (user-accessible) datasets, since they effectively serve as package showcase
  datasets. The following previously exported datasets are therefore downgraded
  to internal datasets, since they are only used in internal functions and
  vignettes to reproduce publication plots: `clim.site.par`, `diffusion.length`
  and `t15.noise`, where the latter replaces the previous combination of the
  `t15` dataset and the internal function `TrenchNoise()`. As a minor change
  alongside that, the T15 trench data, used to create `t15.noise`, is now
  processed directly from the PANGAEA data archive in order to avoid a suggest
  dependency on the `TrenchR` package.

* The input parameter structure to `PlotStackCorrelation()` has been swapped
  such that the function can now directly digest the output of
  `ObtainStackCorrelation()`, e.g. by using a pipe operator.

## New features

* New feature and function `EstimateCI()` that implements a parametric
  bootstrapping procedure to estimate confidence levels for an existing dataset
  of signal, noise, and SNR spectra
  ([#14](https://github.com/EarthSystemDiagnostics/proxysnr/pull/14)).
* `SeparateSignalFromNoise()` now includes the new function parameter
  `measurement.noise` to pass a measurement noise level for correcting the proxy
  noise spectrum, either as a single value (white noise) or as a spectral object
  ([#13](https://github.com/EarthSystemDiagnostics/proxysnr/pull/13)).
* Package documentation now includes the topic `?spec.object`, explaining the
  meaning and expected structure of a `proxysnr` spectral object. This term is
  now consistently used throughout the package and function documentation (see
  also [#12](https://github.com/EarthSystemDiagnostics/proxysnr/pull/12)).
* The functions `LPlot()` and `LLines()` to plot/add a spectrum on a
  double-logarithmic scale are now available as exported functions.
* A new `vignette("proxysnr")` with an introduction to the `proxysnr` package is
  now available.
* A [github.io homepage](https://earthsystemdiagnostics.github.io/proxysnr/) for
  `proxysnr` is now available, created via `pkgdown`. In the process of creating
  that page, the readme, all vignettes, and all function documentation have been
  revised to render nicely on the homepage
  ([#16](https://github.com/EarthSystemDiagnostics/proxysnr/pull/16)).

## Minor updates and changes

* `ObtainArraySpectra()` now accepts the input proxy data both as a list or as a
  data frame.
* `ObtainArraySpectra()`, `SeparateSignalFromNoise()`: the output list of both
  functions now possesses the new attribute `array.par` storing information
  on the original proxy record array (number of proxy records, number of
  observations (e.g. time steps) and (time) resolution), which makes these
  information more easily available, especially for their use in `EstimateCI()`
  ([#15](https://github.com/EarthSystemDiagnostics/proxysnr/pull/15)).
* The `magrittr` pipe operator is now imported into the `proxysnr` namespace,
  since with it several `proxysnr` functions can be executed consecutively;
  e.g. `ObtainArraySpectra() %>% PlotArraySpectra()`, `ObtainArraySpectra() %>%
  SeparateSignalFromNoise()`, `ObtainStackCorrelation() %>%
  PlotStackCorrelation()` (see also
  [#3](https://github.com/EarthSystemDiagnostics/proxysnr/issues/3)).
* `PlotStackCorrelation()`: The color palette parameter has now a default value
  of `NULL`, signalling to set a default color palette internally (taken from the
  ColorBrewer 2.0 collection), which makes use of the function easier and
  quicker.
* `GetIntegratedSNR()`, `ObtainStackCorrelation()`: setting of the integration
  limits has been simplified with the new options to either supply a lower and upper
  frequency axis index, or directly a frequency range.
* `PlotTF()`: It is now possible to only plot a single transfer function type,
  e.g. only a diffusion transfer function or only a time uncertainty transfer
  function ([#5](https://github.com/EarthSystemDiagnostics/proxysnr/pull/5)).
* `CalculateDiffusionTF()`, `CalculateTimeUncertaintyTF()`:
  * the output variables from both functions now include an attribute with
    information on the transfer function simulation run;
  * log-smoothing of the simulated transfer functions can now be switched on via
    the new function parameter `df.log`.
* `PlotArraySpectra()`: Removal of a certain number of spectral estimates on the
  low and/or high-frequency side is no longer hard-coded internally in the
  function but now settable via the new function parameter `remove`, with its
  default setting reproducing the previous default behaviour.
* The internal function `SetPlotPar()` to set default plotting parameters has
  been classified as obsolete and removed.
* Changing the default plotting parameter for the internal plotting functions
  that create the paper plots has been simplified.
* `PublicationSNR()`: Further modified data input arrangement to make function
    work better in conjunction with `WrapSpectralResults()` and to allow
    selection of the data correction version (new function parameter `data`).
* Various small bug fixes.

## Technical updates

* Availability of the suggested package `simproxyage` is now properly being
  checked at the relevant places in the code, and detailed user information how
  to install the package are issued if it is not yet installed.
* The random seed in the vignette on time uncertainty transfer function
  calculation is now correctly set to ensure reproducibe results in every
  dataset case.
* Unit tests have been added, including an automatic GitHub coverage workflow
  (current code coverage > 97%)
  ([#4](https://github.com/EarthSystemDiagnostics/proxysnr/pull/4)).
* The `README` now includes tags for GitHub and Zenodo package versions, the `R
  CMD check` status and the code coverage, and it has been updated to reflect
  the latest status of package availability and documentation.
* Source code now uses 2-space indentation throughout, and all roxygen tags
  now use only the `#'` flag.
* Update to the raw WAIS data processing: the previous automatic download of the
  original published WAIS data has broken due to changes on the data repository
  side. Therefore, to ensure reproducibility, the corresponding data has now
  been manually downloaded and placed into the `data-raw` folder, from where it
  is processed into the package dataset, which exactly reproduces the already
  existing exported dataset in the variable `wais`. For this reason, and to
  preserve the original processing date attribute which reflects the time of the
  corresponding paper publication, the `wais` dataset variable was not
  recreated.
* Vignettes are no longer stored as rendered files under `doc/` since that is
  not supported by `R CMD CHECK`.
* Most utility and paper plotting functions are no longer documented.
* Package license and contributor information has been updated.

# proxysnr 0.2.3

[#257bddb](https://github.com/EarthSystemDiagnostics/proxysnr/commit/257bddba416dc83f6270351dd557152d5e60413c)
* `WrapSpectralResults()`: output element `f.cutoff` with the cutoff frequency
  for constraining the diffusion correction is now found within each output
  element of the raw or corrected spectra, but it has value `NA` for the
  elements where diffusion has not been corrected for.
* Changed data input arrangement of internal function `PublicationSNR()`.

# proxysnr 0.2.2

[#fa9171b](https://github.com/EarthSystemDiagnostics/proxysnr/commit/fa9171bcb0089b6b9a00c59f80487cac79cceddb)
* `DiffusionTF()`: Added new function parameter `window` for analysing only a
  subset of the total simulated record length. Function documentation updated.

# proxysnr 0.2.1

[#17d5ecb](https://github.com/EarthSystemDiagnostics/proxysnr/commit/17d5ecbce68d52e2d7a113f8fbb1f81605470bd5)
* `DiffusionTF()`: Bug fix in the error check on correct input vector length.

# proxysnr 0.2.0 ([#1](https://github.com/EarthSystemDiagnostics/proxysnr/pull/1))

## New function and user-relevant changes

* The function `IntegratedSNR()` has been added, which takes the job of
  calculating the signal-to-noise ratio from the _integrated_ signal and noise
  spectra, so that it delivers the SNR corresponding to varying temporal record
  resolutions (compared to the SNR at a specific frequency obtained from
  `SeparateSpectra()`, which corresponds to a band-pass filtered
  record). `StackCorrelation()` internally now uses the new function for an
  improved code structure.
* Parameter `N` to `StackCorrelation()` no longer indicates the maximum number
  of records but is treated directly as a vector of record numbers; so any call
  which previously used `N = 5` now needs to use `N  = 1 : 5`.
* The rather verbose output of `StackCorrelation()` has been reduced to the
  minimum needed - a list of two elements: the frequency axis and the
  correlation matrix. Output from `StackCorrelation()` which is piped to the
  plotting function `PlotStackCorrelation()` has to be adapted accordingly.
  
## Minor changes
* All code and documentation has been overhauled so that `R CMD check` passes
  with zero errors, warnings or notes.
* The location of built vignettes in the source package for documentation on
  GitHub was changed to reflect recent `devtools` usage.
* The reference to Muench and Laepple (2018) was updated to reflect the final
  publication.

# proxysnr 0.1.0

* Package version released along with the publication Münch, T. and Laepple, T.:
_What climate signal is contained in decadal to centennial scale isotope
variations from Antarctic ice cores?_, Clim. Past, 2018.

