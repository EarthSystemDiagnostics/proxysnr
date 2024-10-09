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
