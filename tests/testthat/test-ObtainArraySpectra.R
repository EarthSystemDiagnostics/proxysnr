test_that("obtaining the array spectra works", {

  # test error checking

  m <- "`cores` must be a list or a data frame."
  expect_error(ObtainArraySpectra(1), m, fixed = TRUE)
  expect_error(ObtainArraySpectra(matrix(1 : 6, 2, 3)), m, fixed = TRUE)

  m <- "Only one proxy record (`cores` is of length 1)."
  expect_error(ObtainArraySpectra(list(a = 1 : 10)), m, fixed = TRUE)

  m <- "Elements of `cores` must all have the same length."
  expect_error(ObtainArraySpectra(list(a = 1, b = 1 : 10, c = c(3, 8, 5))),
                                  m, fixed = TRUE)
  
  # test with data

  cores <- list(a = rnorm(100), b = rnorm(100), c = rnorm(100))

  cores.ts <- lapply(cores, stats::ts)
  single <- lapply(cores.ts, SpecMTM)
  mean <- MeanSpectrum(single)
  stack <- SpecMTM(stats::ts(rowMeans(simplify2array(cores))))

  expected <- list(single = single, mean = mean, stack = stack)
  attr(expected, "array.par") <- c(nc = 3, nt = 100, res = 1)
  actual <- ObtainArraySpectra(cores)

  expect_equal(actual, expected)
  expect_equal(names(actual$single), names(cores))

  # different frequency axis

  cores.ts <- lapply(cores, stats::ts, deltat = 5)
  single <- lapply(cores.ts, SpecMTM)
  mean <- MeanSpectrum(single)
  stack <- SpecMTM(stats::ts(rowMeans(simplify2array(cores)), deltat = 5))

  expected <- list(single = single, mean = mean, stack = stack)
  attr(expected, "array.par") <- c(nc = 3, nt = 100, res = 5)
  actual <- ObtainArraySpectra(cores, res = 5)

  expect_equal(actual, expected)

  # effective number of records is set

  neff <- 1.5
  expected <- list(single = single, mean = mean, stack = stack)
  attr(expected, "array.par") <- c(nc = neff, nt = 100, res = 5)
  actual <- ObtainArraySpectra(cores, res = 5, neff = neff)

  expect_equal(actual, expected)

  # test deprecated function name

  expect_warning(actual <- ArraySpectra(cores, res = 5, neff = neff))
  expect_equal(actual, expected)

})
