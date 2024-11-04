test_that("SeparateSignalFromNoise error checks work", {

  # --- test main input checks -------------------------------------------------

  m <- "`spectra` must be a list."
  expect_error(SeparateSignalFromNoise(1), m, fixed = TRUE)

  m <- "`spectra` must have elements `mean` and `stack`."
  expect_error(SeparateSignalFromNoise(list(foo = 1)), m, fixed = TRUE)
  expect_error(SeparateSignalFromNoise(list(mean = 1 : 10)), m, fixed = TRUE)
  expect_error(SeparateSignalFromNoise(list(stack = 1 : 10)), m, fixed = TRUE)

  mean <- stack <- list(freq = 1 : 10, spec = 1 : 10)

  m <- "`spectra$mean` must be a list with elements `freq` and `spec` of equal length."
  expect_error(SeparateSignalFromNoise(list(mean = 1, stack = stack)),
               m, fixed = TRUE)
  m <- "`spectra$stack` must be a list with elements `freq` and `spec` of equal length."
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = 1)),
               m, fixed = TRUE)

  m <- "`mean` and `stack` must have the same number of spectral estimates."
  short_stack <- list(freq = 1 : 4, spec = 1 : 4)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = short_stack)),
               m, fixed = TRUE)

  m <- "Frequency axes of `mean` and `stack` do not match."
  wrong_freq <- list(freq = 2 : 11, spec = 1 : 10)
  expect_error(SeparateSignalFromNoise(list(mean = wrong_freq, stack = stack)),
               m, fixed = TRUE)

  m <- "Supply (effective) number of records."
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack)),
               m, fixed = TRUE)

  # --- test input checks related to measurement noise -------------------------

  m <- paste("`measurement.noise` must be a single value or a spectral object.")
  measurement.noise <- rnorm(3)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       measurement.noise = measurement.noise),
               m, fixed = TRUE)

  # should also work with length-1 array
  measurement.noise <- matrix(0.1, 1, 1)
  expect_no_warning(SeparateSignalFromNoise(
    list(mean = mean, stack = stack, N = 1),
    measurement.noise = measurement.noise))

  m <- paste("`measurement.noise` must be a list with elements",
             "`freq` and `spec` of equal length.")
  measurement.noise <- list(a = 1)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       measurement.noise = measurement.noise),
               m, fixed = TRUE)

  m <- paste("No sufficient frequency axis overlap between proxy data",
             "and measurement noise spectrum.")
  measurement.noise <- list(freq = 0.5, spec = 0.1)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       measurement.noise = measurement.noise),
               m, fixed = TRUE)

  # --- test input checks related to transfer functions ------------------------

  m <- paste("`diffusion` must be a list with elements",
             "`freq` and `spec` of equal length.")
  tf <- list(a = 1)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       diffusion = tf), m, fixed = TRUE)
  tf <- list(freq = 1, spec = 1 : 5)
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       diffusion = tf), m, fixed = TRUE)
  m <- paste("`time.uncertainty` must be a list with elements",
             "`freq` and `spec` of equal length.")
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       time.uncertainty = tf), m, fixed = TRUE)

  tf <- list(freq = 1 : 5, spec = 1 : 5)
  m <- paste("No sufficient frequency axis overlap between proxy data",
             "and diffusion transfer function.")
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       diffusion = tf), m, fixed = TRUE)
  m <- paste("No sufficient frequency axis overlap between proxy data",
             "and time uncertainty transfer function.")
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       time.uncertainty = tf), m, fixed = TRUE)

})

test_that("SeparateSignalFromNoise calculations work", {

  # --- test default behaviour -------------------------------------------------

  signal <- list(freq = 1 : 10, spec = 1 : 10)
  noise <- list(freq = 1 : 10, spec = rep(5, 10))
  snr <- list(freq = 1 : 10, spec = signal$spec / noise$spec)
  class(signal) <- class(noise) <- class(snr) <- "spec"

  n <- 5
  mean <- list(freq = signal$freq, spec = signal$spec + noise$spec)
  stack <- list(freq = signal$freq, spec = signal$spec + noise$spec / n)

  actual <- SeparateSignalFromNoise(spectra = list(mean = mean, stack = stack),
                                    neff = n)
  expected <- list(N = n, signal = signal, noise = noise, snr = snr)

  expect_equal(actual, expected)

  # --- test deprecated function name ------------------------------------------

  expect_warning(actual <- SeparateSpectra(spectra = list(mean = mean, stack = stack),
                                           neff = n))
  expect_equal(actual, expected)

  # --- test including N and transfer function input ---------------------------

  diff <- list(freq = signal$freq, spec = rep(1 / 1.5, length(signal$freq)))
  tunc <- list(freq = signal$freq, spec = rep(1 / 1.1, length(signal$freq)))

  mean <- list(freq = signal$freq,
               spec = diff$spec * (signal$spec + noise$spec))
  stack <- list(freq = signal$freq,
                spec = diff$spec * (tunc$spec * signal$spec + noise$spec / n))

  actual <- SeparateSignalFromNoise(
    spectra = list(mean = mean, stack = stack, N = n),
    diffusion = diff,
    time.uncertainty = tunc)
  expected <- list(N = n, signal = signal, noise = noise, snr = snr)

  expect_equal(actual, expected)

  # --- test interpolation of transfer functions -------------------------------

  data <- ObtainArraySpectra(dml$dml2)

  dtf <- diffusion.tf$dml2
  ttf <- time.uncertainty.tf$dml2
  expected <- data %>%
    SeparateSignalFromNoise(diffusion = dtf, time.uncertainty = ttf)

  # artificially extend transfer functions outside range of data
  # <- extension part should be exactly removed by interpolation routine
  dtf.ext <- list(
    freq = c(0, 0.0005, 0.0006, 0.0009, dtf$freq, 0.5, 0.6, 0.7, 0.9, 1, 10),
    spec = c(runif(4, 0.1, 1), dtf$spec, runif(6, 0.1, 1))
  )
  ttf.ext <- list(
    freq = c(0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, ttf$freq, 1, 10),
    spec = c(runif(6, 0.1, 1), ttf$spec, runif(2, 0.1, 1))
  )

  actual <- data %>%
    SeparateSignalFromNoise(diffusion = dtf.ext, time.uncertainty = ttf.ext)

  expect_equal(actual, expected)

  # test applying transfer function calculated on semi-annual resolution on
  # annual resolution package data
  tf <- CalculateDiffusionTF(nt = length(1994 : 1000) / 0.5, nc = 3, ns = 3,
                             sigma = proxysnr:::diffusion.length$dml2[, -1],
                             res = 0.5)

  expect_no_error(SeparateSignalFromNoise(data, diffusion = tf))

  # --- test measurement noise input -------------------------------------------

  # PSD measurement noise level of a white spectrum
  measurement_noise_psd <- 0.25

  # white measurement noise given as a single value

  signal <- list(freq = 1 : 10, spec = 1 : 10)
  noise <- list(freq = 1 : 10, spec = rep(5, 10))
  snr <- list(freq = 1 : 10, spec = signal$spec / noise$spec)
  class(signal) <- class(noise) <- class(snr) <- "spec"

  n <- 5
  mean <- list(freq = signal$freq,
               spec = signal$spec + noise$spec + measurement_noise_psd)
  stack <- list(freq = signal$freq,
                spec = signal$spec + noise$spec / n + measurement_noise_psd / n)

  # measurement noise variance
  measurement_noise_var <- 2 * max(signal$freq) * measurement_noise_psd
  actual <- SeparateSignalFromNoise(spectra = list(mean = mean, stack = stack),
                                    neff = n,
                                    measurement.noise = measurement_noise_var)
  expected <- list(N = n, signal = signal, noise = noise, snr = snr)

  expect_equal(actual, expected)

  # white measurement noise given as a spectral object

  measurement_noise <- list(freq = c(0.5, 3, 8, 15),
                            spec = rep(measurement_noise_psd, 4))

  actual <- SeparateSignalFromNoise(spectra = list(mean = mean, stack = stack),
                                    neff = n, measurement.noise = measurement_noise)

  expect_equal(actual, expected)

})
