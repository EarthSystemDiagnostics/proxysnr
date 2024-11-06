signal.par <- list(alpha = 0.1, beta = 1)
noise.par  <- list(alpha = 0.05, beta = 0.1)

nc <- 5
nt <- 100

nmc <- 4

test_that("simulation of proxy core array works", {

  sim <- simCoreArray(signal.par, noise.par, nc = 5, nt = 100)

  expect_type(sim, "list")
  expect_length(sim, nc)
  expect_equal(lengths(sim), rep(nt, nc))

})

test_that("simulated signal and noise spectra are valid", {
  
  sim <- simSignalAndNoise(signal.par, noise.par, nc = nc, nt = nt, res = 1)

  expect_type(sim, "list")
  expect_length(sim, 4)
  expect_named(sim, c("N", "signal", "noise", "snr"))

  expect_true(is.spectrum(sim$signal))
  expect_true(is.spectrum(sim$noise))
  expect_true(is.spectrum(sim$snr))

})

test_that("running surrogate signal and noise spectra works", {

  sim <- runSurrogates(signal.par, noise.par, nc = nc, nt = nt,
                       res = 1, nmc = nmc)

  expect_type(sim, "list")
  expect_length(sim, nmc)
  expect_equal(lengths(sim), rep(4, nmc))

})

test_that("extracting relative quantiles from realizations works", {

  foo <- function(probs = c(0.1, 0.9), res = 1) {
    runSurrogates(signal.par, noise.par, nc = nc, nt = nt,
                      res = res, nmc = nmc) %>%
      extractQuantiles(probs)
  }

  m <- "`probs` must be a length-2 vector."
  expect_error(foo(probs = 0.1), m, fixed = TRUE)
  expect_error(foo(probs = c(0.1, 0.5, 0.9)), m, fixed = TRUE)

  m <- "`probs[2]` must be > `probs[1]`."
  expect_error(foo(probs = c(0.9, 0.1)), m, fixed = TRUE)
  expect_error(foo(probs = c(0.1, 0.1)), m, fixed = TRUE)

  ci <- foo()

  # expected freq axis (valid for deltat = 1, else scale by 1 / deltat)
  f <- seq(1 / nt, 0.5, by = 1 / nt)
  nf <- length(f)

  expect_type(ci, "list")
  expect_length(ci, 3)
  expect_named(ci, c("signal", "noise", "snr"))

  expect_named(ci$signal, c("freq", "lower", "upper"))
  expect_equal(lengths(ci$signal, use.names = FALSE), rep(nf, 3))
  expect_equal(ci$signal$freq, f)

  expect_named(ci$noise, c("freq", "lower", "upper"))
  expect_equal(lengths(ci$noise, use.names = FALSE), rep(nf, 3))
  expect_equal(ci$noise$freq, f)

  expect_named(ci$snr, c("freq", "lower", "upper"))
  expect_equal(lengths(ci$snr, use.names = FALSE), rep(nf, 3))
  expect_equal(ci$snr$freq, f)

  expect_true(all(ci$signal$lower < 1))
  expect_true(all(ci$signal$upper > 1))

  expect_true(all(ci$noise$lower < 1))
  expect_true(all(ci$noise$upper > 1))

  expect_true(all(ci$snr$lower < 1))
  expect_true(all(ci$snr$upper > 1))

  # check with non-unity deltat

  deltat <- 1.7
  ci <- foo(res = deltat)
  f <- seq(1 / nt / deltat, 0.5 / deltat, by = 1 / nt / deltat)
  nf <- length(f)

  expect_equal(lengths(ci$signal, use.names = FALSE), rep(nf, 3))
  expect_equal(ci$signal$freq, f)

})

test_that("simulation produces correct frequency axis", {

  # with default deltat (res) setting
  spectra <- ObtainArraySpectra(dml$dml2, df.log = 0.15) %>%
    SeparateSignalFromNoise()

  sim <- runSimulation(spectra, f.end = 0.1, nmc = 1)

  expect_equal(spectra$signal$freq, sim[[1]]$signal$freq)

  # with some custom deltat (res) setting
  res <- 3.89
  spectra <- ObtainArraySpectra(dml$dml2, df.log = 0.15, res = res) %>%
    SeparateSignalFromNoise()

  sim <- runSimulation(spectra, f.end = 0.1, nmc = 1)

  expect_equal(spectra$signal$freq, sim[[1]]$signal$freq)

})

test_that("CI error checks work", {

  m <- "`spectra` must be a list."

  expect_error(EstimateCI(1), m, fixed = TRUE)
  expect_error(EstimateCI(matrix()), m, fixed = TRUE)

  m <- "`spectra` must have elements `N`, `signal`, `noise`, and `snr`."

  expect_error(EstimateCI(list()), m, fixed = TRUE)
  expect_error(EstimateCI(list(N = 1)), m, fixed = TRUE)
  expect_error(EstimateCI(list(N = 1, signal = list())), m, fixed = TRUE)
  expect_error(EstimateCI(list(signal = list(), noise = list())),
               m, fixed = TRUE)
  expect_error(EstimateCI(list(N = 1, signal = list(), snr = list())),
               m, fixed = TRUE)

  s1 <- list(freq = 1, spec = 1)
  s2 <- list(freq = 3 : 12, spec = rnorm(10))

  m <- paste("`spectra$signal` must be a list with elements",
             "`freq` and `spec` of equal length.")

  expect_error(EstimateCI(
    list(N = 1, signal = list(), noise = list(), snr = list())),
    m, fixed = TRUE)

  m <- paste("`spectra$snr` must be a list with elements",
             "`freq` and `spec` of equal length.")

  expect_error(EstimateCI(
    list(N = 1, signal = s1, noise = s1, snr = list())),
    m, fixed = TRUE)

  m <- "`signal` and `noise` must have the same number of spectral estimates."

  expect_error(EstimateCI(
    list(N = 1, signal = s1, noise = s2, snr = s1)),
    m, fixed = TRUE)

  m <- "`signal` and `snr` must have the same number of spectral estimates."

  expect_error(EstimateCI(
    list(N = 1, signal = s1, noise = s1, snr = s2)),
    m, fixed = TRUE)

  s1 <- list(freq = 1 : 10, spec = rnorm(10))

  m <- "Frequency axes of `signal` and `noise` do not match."

  expect_error(EstimateCI(
    list(N = 1, signal = s1, noise = s2, snr = s2)),
    m, fixed = TRUE)

  m <- "Frequency axes of `signal` and `snr` do not match."

  expect_error(EstimateCI(
    list(N = 1, signal = s1, noise = s1, snr = s2)),
    m, fixed = TRUE)

  m <- paste("`spectra$N` must be a single integer > 1 supplying the number",
             "of records underlying the analyses in `spectra`.")

  spectra <- list(signal = s1, noise = s1, snr = s1)

  spectra$N <- 1 : 3
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- -5
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- 1
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- list(N = 1)
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- data.frame(N = 1)
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- matrix(data = 1, ncol = 2)
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- NA
  expect_error(EstimateCI(spectra), m, fixed = TRUE)
  spectra$N <- Inf
  expect_error(EstimateCI(spectra), m, fixed = TRUE)

})

test_that("CI estimation works", {

  # test using DML2 dataset from Muench and Laepple (2018)

  spectra <- ObtainArraySpectra(dml$dml2, df.log = 0.15) %>%
    SeparateSignalFromNoise(diffusion = diffusion.tf$dml2)

  spectra.with.ci <- EstimateCI(spectra, f.end = 0.1, nmc = 3,
                                df.log = NULL, ci.df.log = NULL)

  expect_type(spectra.with.ci, "list")
  expect_length(spectra.with.ci, 4)
  expect_named(spectra.with.ci, c("N", "signal", "noise", "snr"))

  nms <- c("freq", "spec", "lim.1", "lim.2")

  expect_type(spectra.with.ci$signal, "list")
  expect_true(all(utils::hasName(spectra.with.ci$signal, nms)))

  expect_type(spectra.with.ci$noise, "list")
  expect_true(all(utils::hasName(spectra.with.ci$noise, nms)))

  expect_type(spectra.with.ci$snr, "list")
  expect_true(all(utils::hasName(spectra.with.ci$snr, nms)))

  # test same dataset with log-smoothing switched on

  expect_no_error(
    EstimateCI(spectra, f.end = 0.1, nmc = 3, df.log = 0.1, ci.df.log = 0.05))

  # test the combined DML1+DML2 dataset, which has a merged frequency axis, so
  # that interpolation of the simulated CI should be needed

  w <- paste0("Input and simulated frequency axes differ;\n",
              "confidence intervals are interpolated possibly ",
              "producing NA values.")

  spectra <- WrapSpectralResults(
    dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
    diffusion = diffusion.tf,
    time.uncertainty = time.uncertainty.tf,
    df.log = c(0.15, 0.15, 0.1)) %>%
    proxysnr:::PublicationSNR(data = "raw") %>%
    .$dml

  # supply number of records manually: use 3
  spectra$N <- 3

  expect_warning(
    spectra.with.ci <- EstimateCI(spectra, f.end = 0.1,
                                  df.log = 0.15, ci.df.log = 0.05),
    w)

  expect_type(spectra.with.ci, "list")
  expect_length(spectra.with.ci, 5)
  expect_named(spectra.with.ci, c("signal", "noise", "snr", "f.cutoff", "N"))
  expect_equal(spectra.with.ci$N, spectra$N)

  nms <- c("freq", "spec", "lim.1", "lim.2")

  expect_type(spectra.with.ci$signal, "list")
  expect_true(all(utils::hasName(spectra.with.ci$signal, nms)))

  expect_type(spectra.with.ci$noise, "list")
  expect_true(all(utils::hasName(spectra.with.ci$noise, nms)))

  expect_type(spectra.with.ci$snr, "list")
  expect_true(all(utils::hasName(spectra.with.ci$snr, nms)))


})
