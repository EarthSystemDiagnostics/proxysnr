#' Time uncertainty transfer function
#' 
#' This function implements an empirical Monte Carlo approach to estimate the
#' spectral transfer function for the effect of time uncertainty on the spatial
#' average of a common proxy signal recorded by a given core array. It requires
#' the suggested package \code{simproxyage}; if the package is not installed,
#' you will receive installation instructions upon calling this function.
#'
#' The approach is described in detail in Münch and Laepple (2018). In brief,
#' \code{nc} identical surrogate time series are created and age perturbed and
#' the average of these time series is calculated. The process is repeated
#' \code{ns} times. For each of the \code{ns} realisations, spectra of the
#' average age-perturbed and original time series are calculated yielding the
#' spectral transfer function.
#'
#' The modelling of the time uncertainty follows the approach presented in
#' Comboul et al. (2014) which is implemented in the package
#' \code{simproxyage}. The package can be installed from GitHub via\cr
#' `remotes::install_github("EarthSystemDiagnostics/simproxyage")`\cr
#' or you can download the package repository from\cr
#' <https://doi.org/10.5281/zenodo.2025204>\cr
#' and install the package via `devtools::install()`.
#'
#' The spectral estimates are calculated using Thomson’s multitaper method with
#' three windows with linear detrending before analysis.
#'
#' Handling of age control points (see package \code{simproxyage} for more
#' details):
#' 
#' The function allows one to model the age uncertainty with several age
#' control points (ACP) where the age uncertainty evolves according to the
#' model of Comboul et al. (2014) between the age control points but is forced
#' back to zero time uncertainty at the ACPs following a Brownian Bridge
#' process.
#'
#' Per default, the start age (i.e. the youngest age at the core top) is
#' assumed to be the first ACP, thus, if not set explicitly in the vector
#' \code{acp}, the youngest age is added as its first element. You can specify
#' an arbitrary number of additional ACPs. Between each pair of ACPs (starting
#' at the core top), constrained age perturbation realisations are performed
#' following the Comboul et al. (2014) model and the Brownian Bridge
#' concept. If the last ACP equals the oldest age in \code{t} (i.e. the last
#' proxy data point), the last core segment also follows a Brownian Bridge
#' process. Alternatively, if the last ACP is younger than the final age,
#' \code{NA} is added as the last element to the vector \code{acp} which
#' results in an unconstrained age perturbation process for the last core
#' segment.
#'
#' ACPs lying outside the time interval defined by \code{t} will be removed
#' from the vector \code{acp} with a warning.
#'
#' @param t numeric vector of integer values providing a reference chronology
#'   for the age perturbations (starting with the youngest age).
#' @param acp numeric vector of age control points where the age uncertainty is
#'   assumed to be zero. Per default, a two-element vector where the first
#'   element is the start age (\code{t[1]}) and the second element \code{NA}
#'   resulting in an unconstrained age perturbation process; see details for
#'   handling more age control points.
#' @param nt the length of the records (i.e. the number of data points) to
#'   simulate; per default set to \code{length(t)}.
#' @param nc the number of cores in the modelled core array.
#' @param ns the number of Monte Carlo simulations for estimating the transfer
#'   function.
#' @param model name string of the random process to use for realising the age
#'   perturbations; must be either "poisson" (the default) or "bernoulli"; see
#'   Comboul et al. (2014) for details on the two models.
#' @param rate numeric vector of probability rate(s) that an age band is
#'   perturbed; you can specify a vector of two rates where the first entry is
#'   the probability for a missing band and the second entry the probability for
#'   a double-counting of a band. If only a single value is specified (per
#'   default 0.05), symmetric perturbations are assumed.
#' @param resize the resizing option in case of shorter/longer than original
#'   time axes: 0 = do not resize, -1 = resize to shortest realisation, 1 =
#'   resize to longest realisation (default).
#' @param surrogate.fun the random number generator to use for creating the
#'   noise time series; per default, Gaussian white noise is created using the
#'   base \code{R} function \code{\link{rnorm}}.
#' @param fun.par an otional list of additional parameters which are passed to
#'   \code{surrogate.fun}.
#' @param pad Strong age perturbations may result in time series that do not
#'   reach the final age of the reference chronology -- shall the resulting
#'   \code{NA} values (in case of \code{resize = 1}) be set to \code{0} for the
#'   spectral estimation? Defaults to \code{TRUE}. This should not influence the
#'   spectral estimation provided the number of \code{NA} values is small
#'   compared to the total length of the time series.
#' @param ... additional parameters which are passed to the spectral estimation
#'   function \code{\link{SpecMTM}}.
#'
#' @return a list of the components \code{input}, \code{stack} and
#'   \code{ratio} which are objects of class \code{"spec"} providing averages
#'   over the \code{ns} simulations of:
#'   \describe{
#'   \item{\code{input}:}{the original spectrum of the surrogate data;}
#'   \item{\code{stack}:}{the spectrum of the spatial average of the
#'         age-perturbed records;}
#'   \item{\code{ratio}:}{their ratio (perturbed/unperturbed), i.e. the transfer
#'         function.}
#' }
#'
#' @author Thomas Münch
#'
#' @references
#'
#' Comboul, M., Emile-Geay, J., Evans, M. N., Mirnateghi, N., Cobb, K. M.
#' and Thompson, D. M.: A probabilistic model of chronological errors in
#' layer-counted climate proxies: applications to annually banded coral
#' archives, Clim. Past, 10(2), 825-841,
#' https://doi.org/10.5194/cp-10-825-2014, 2014.
#'
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, https://doi.org/10.5194/cp-14-2053-2018, 2018.
#'
#' @export
#'
CalculateTimeUncertaintyTF <- function(t = 100 : 1, acp = c(t[1], NA),
                                       nt = length(t), nc = 1, ns = 100,
                                       model = "poisson", rate = 0.05,
                                       resize = 1, surrogate.fun = stats::rnorm,
                                       fun.par = NULL, pad = TRUE, ...) {

  # check if package simproxyage is available
  has.simproxyage <- check.simproxyage(stop.on.false = TRUE)

  # Monte Carlo simulation of age uncertainty for a core array
  arg <- c(list(
    t = t,
    acp = acp,
    nt = nt,
    nc = nc,
    ns = ns,
    model = model,
    rate = rate,
    resize = resize,
    surrogate.fun = surrogate.fun),
    fun.par)

  run <- do.call(simproxyage::MonteCarloArray, args = arg)
  
  
  # pad potential NA's at the end of the age-perturbed time series with zero
  stacks <- run$stacks
  if (pad) {
    stacks[which(is.na(run$stacks))] <- 0
  }


  # calculate spectra of the 'ns' input and age-perturbed time series

  stacks.spec <- lapply(seq_len(ncol(stacks)), function(i) {
    SpecMTM(stats::ts(stacks[, i], deltat = 1), ...)})

  input.spec <- lapply(seq_len(ncol(run$input)), function(i) {
    SpecMTM(stats::ts(run$input[, i], deltat = 1), ...)})


  # calculate the average over all 'ns' simulations

  stacks.spec.mean <- MeanSpectrum(stacks.spec)
  input.spec.mean  <- MeanSpectrum(input.spec)


  # return average input and age-perturbed spectra and the corresponding
  # ratio (transfer function) as a list

  res <- list()
  res$input <- input.spec.mean
  res$stack <- stacks.spec.mean
  res$ratio$freq <- stacks.spec.mean$freq
  res$ratio$spec <- stacks.spec.mean$spec / input.spec.mean$spec

  class(res$ratio) <- "spec"

  return(res)

}

#' Check if `simproxyage` package is installed
#'
#' @param stop.on.false logical whether to stop when package is not installed;
#'   defaults to `FALSE` which only prints the message with the installation
#'   instructions.
#' @return a logical signalling whether the package is installed.
#'
#' @author Thomas Münch
#' @noRd
#'
check.simproxyage <- function(stop.on.false = FALSE) {

  has.simproxyage <- requireNamespace("simproxyage", quietly = TRUE)

  if (has.simproxyage) {

    message("\nChecking simproxyage availability... ok.")

  } else {

    message("\nChecking simproxyage availability...")

    msg <- paste0(
      "Package not found.\n",
      "Install via ",
      "`remotes::install_github(\"EarthSystemDiagnostics/simproxyage\")`\n",
      "or download package repository from ",
      "<https://doi.org/10.5281/zenodo.2025204>\n",
      "and install via `devtools::install()`.\n"
    )

    if (stop.on.false) {
      cat("\n")
      msg <- paste0(msg, "\n")
      stop(msg, call. = FALSE)
    } else {
      msg <- paste0("\n", msg)
      message(msg)
    }
  }

  invisible(has.simproxyage)

}
