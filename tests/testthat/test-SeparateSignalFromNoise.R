test_that("SeparateSignalFromNoise error checks work", {

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
  m <- paste("Length of diffusion correction does not match",
             "length of spectral estimates.")
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       diffusion = tf), m, fixed = TRUE)
  m <- paste("Length of time uncertainty correction does not match",
             "length of spectral estimates.")
  expect_error(SeparateSignalFromNoise(list(mean = mean, stack = stack, N = 1),
                                       time.uncertainty = tf), m, fixed = TRUE)

})

test_that("SeparateSignalFromNoise calculations work", {

  signal <- list(freq = 1 : 10, spec = 1 : 10)
  noise <- list(freq = 1 : 10, spec = rep(5, 10))
  snr <- list(freq = 1 : 10, spec = signal$spec / noise$spec)
  class(signal) <- class(noise) <- class(snr) <- "spec"

  n <- 5
  mean <- list(freq = signal$freq, spec = signal$spec + noise$spec)
  stack <- list(freq = signal$freq, spec = signal$spec + noise$spec / n)

  actual <- SeparateSignalFromNoise(spectra = list(mean = mean, stack = stack),
                                    neff = n)
  expected <- list(signal = signal, noise = noise, snr = snr)

  expect_equal(actual, expected)

  # test deprecated function name

  expect_warning(actual <- SeparateSpectra(spectra = list(mean = mean, stack = stack),
                                           neff = n))
  expect_equal(actual, expected)

  # with diffusion and time uncertainty correction and N included in input

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
  expected <- list(signal = signal, noise = noise, snr = snr)

  expect_equal(actual, expected)

})
