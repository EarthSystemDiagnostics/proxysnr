# proxysnr 0.2.1

* `DiffusionTF()`: Bug fix in the error check on correct input vector length.

# proxysnr 0.2.0

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

