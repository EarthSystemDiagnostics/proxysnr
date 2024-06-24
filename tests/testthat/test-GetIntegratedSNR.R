test_that("GetIntegratedSNR error checks work", {

  m <- "`input` must be a list."
  expect_error(GetIntegratedSNR(1), m, fixed = TRUE)

  m <- "`input` must have elements `signal` and `noise`."
  expect_error(GetIntegratedSNR(list(foo = 1)), m, fixed = TRUE)
  expect_error(GetIntegratedSNR(list(mean = 1 : 10)), m, fixed = TRUE)
  expect_error(GetIntegratedSNR(list(stack = 1 : 10)), m, fixed = TRUE)

  signal <- noise <- list(freq = 1 : 10, spec = 1 : 10)

  m <- "`input$signal` must be a list with elements `freq` and `spec` of equal length."
  expect_error(GetIntegratedSNR(list(signal = 1, noise = noise)),
               m, fixed = TRUE)
  m <- "`input$noise` must be a list with elements `freq` and `spec` of equal length."
  expect_error(GetIntegratedSNR(list(signal = signal, noise = 1)),
               m, fixed = TRUE)

  m <- "`signal` and `noise` must have the same number of spectral estimates."
  short_signal <- list(freq = 1 : 4, spec = 1 : 4)
  expect_error(GetIntegratedSNR(list(signal = short_signal, noise = noise)),
               m, fixed = TRUE)

  m <- "Frequency axes of `signal` and `noise` do not match."
  wrong_freq <- list(freq = 2 : 11, spec = 1 : 10)
  expect_error(GetIntegratedSNR(list(signal = signal, noise = wrong_freq)),
               m, fixed = TRUE)

  m <- "`N` must have length 1."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), N = 1 : 5),
    m, fixed = TRUE)
  m <- "`N` must be >= 1."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), N = -1),
    m, fixed = TRUE)

  m <- "`limits` must have length 2."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), limits = 1 : 5),
    m, fixed = TRUE)
  m <- "`limits[2]` must be > `limits[1]`."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), limits = c(6, 4)),
    m, fixed = TRUE)

  m <- "`f1` must have length 1."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f1 = 1 : 6),
    m, fixed = TRUE)
  m <- "`f2` must have length 1."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f2 = 1 : 6),
    m, fixed = TRUE)

  m <- "`f1` must be >= 1."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f1 = -1),
    m, fixed = TRUE)

  m <- "`f2` must be set to \"max\" or to an index."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f2 = "foo"),
    m, fixed = TRUE)

  m <- "`f2` larger than length of data."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f2 = 13),
    m, fixed = TRUE)

  m <- "`f2` must be > `f1`."
  expect_error(
    GetIntegratedSNR(list(signal = signal, noise = noise), f1 = 6, f2 = 4),
    m, fixed = TRUE)

})

test_that("integration for snr calculation works", {

  signal <- noise <- list(freq = 1 : 10, spec = 1 : 10)
  input <- list(signal = signal, noise = noise)

  expected <- list(freq = 2 : 10, spec = rep(1, 9))
  class(expected) <- "spec"
  actual <- GetIntegratedSNR(input)

  expect_equal(actual, expected)

  # test deprecated function name
  expect_warning(actual <- IntegratedSNR(input))
  expect_equal(actual, expected)

  expected <- list(freq = 2 : 10, spec = rep(4, 9))
  class(expected) <- "spec"
  actual <- GetIntegratedSNR(input, N = 4)

  expect_equal(actual, expected)

  expected <- list(freq = 4 : 6, spec = rep(1, 3))
  class(expected) <- "spec"
  actual <- GetIntegratedSNR(input, f1 = 4, f2 = 6)

  expect_equal(actual, expected)

  expected <- list(freq = 1 : 7, spec = rep(1, 7))
  class(expected) <- "spec"
  actual <- GetIntegratedSNR(input, limits = c(1, 7))

  expect_equal(actual, expected)

  # test deprecated function name
  expect_warning(actual <- IntegratedSNR(input, freq.cut.lower = 1,
                                         freq.cut.upper = 7))
  expect_equal(actual, expected)

})
