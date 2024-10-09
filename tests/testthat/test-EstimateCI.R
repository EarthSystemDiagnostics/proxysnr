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
  expect_length(sim, 3)
  expect_named(sim, c("signal", "noise", "snr"))

  expect_true(is.spectrum(sim$signal))
  expect_true(is.spectrum(sim$noise))
  expect_true(is.spectrum(sim$snr))

})

test_that("running surrogate signal and noise spectra works", {

  sim <- runSurrogates(signal.par, noise.par, nc = nc, nt = nt,
                       res = 1, nmc = nmc)

  expect_type(sim, "list")
  expect_length(sim, nmc)
  expect_equal(lengths(sim), rep(3, nmc))

})

test_that("extracting relative quantiles from realizations works", {

  foo <- function(probs = c(0.1, 0.9)) {
    runSurrogates(signal.par, noise.par, nc = nc, nt = nt,
                      res = 1, nmc = nmc) %>%
      extractQuantiles(probs)
  }

  m <- "`probs` must be a length-2 vector."
  expect_error(foo(probs = 0.1), m, fixed = TRUE)
  expect_error(foo(probs = c(0.1, 0.5, 0.9)), m, fixed = TRUE)

  m <- "`probs[2]` must be > `probs[1]`."
  expect_error(foo(probs = c(0.9, 0.1)), m, fixed = TRUE)
  expect_error(foo(probs = c(0.1, 0.1)), m, fixed = TRUE)

  ci <- foo()

  expect_type(ci, "list")
  expect_length(ci, 3)
  expect_named(ci, c("signal", "noise", "snr"))

  expect_named(ci$signal, c("lower", "upper"))
  expect_equal(lengths(ci$signal, use.names = FALSE), rep(nt / 2, 2))

  expect_named(ci$noise, c("lower", "upper"))
  expect_equal(lengths(ci$noise, use.names = FALSE), rep(nt / 2, 2))

  expect_named(ci$snr, c("lower", "upper"))
  expect_equal(lengths(ci$snr, use.names = FALSE), rep(nt / 2, 2))

  expect_true(all(ci$signal$lower < 1))
  expect_true(all(ci$signal$upper > 1))

  expect_true(all(ci$noise$lower < 1))
  expect_true(all(ci$noise$upper > 1))

  expect_true(all(ci$snr$lower < 1))
  expect_true(all(ci$snr$upper > 1))

})

test_that("CI estimation works", {

  spectra <- ObtainArraySpectra(dml$dml2, df.log = 0.15) %>%
    SeparateSignalFromNoise(diffusion = diffusion.tf$dml2)

  spectra.with.ci <- EstimateCI(spectra, f.end = 0.1, nc = 3, res = 1,
                                nmc = 3, df.log = NULL, ci.df.log = NULL)

  expect_type(spectra.with.ci, "list")
  expect_length(spectra.with.ci, 3)
  expect_named(spectra.with.ci, c("signal", "noise", "snr"))

  nms <- c("freq", "spec", "lim.1", "lim.2")

  expect_type(spectra.with.ci$signal, "list")
  expect_true(all(utils::hasName(spectra.with.ci$signal, nms)))

  expect_type(spectra.with.ci$noise, "list")
  expect_true(all(utils::hasName(spectra.with.ci$noise, nms)))

  expect_type(spectra.with.ci$snr, "list")
  expect_true(all(utils::hasName(spectra.with.ci$snr, nms)))

  # test with log-smoothing switched on
  expect_no_error(
    EstimateCI(spectra, f.end = 0.1, nc = 3, res = 1,
               nmc = 3, df.log = 0.1, ci.df.log = 0.05))

})
