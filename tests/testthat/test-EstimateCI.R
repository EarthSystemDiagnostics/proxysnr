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
