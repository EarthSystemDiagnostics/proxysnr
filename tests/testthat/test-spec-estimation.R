test_that("spectral estimation works", {

  m <- "missing data"
  expect_error(SpecMTM(stats::ts(c(rnorm(41), rep(NA, 9), rnorm(50)))), m)

  # only simple testing here on correct output structure

  data <- rnorm(1000)

  actual <- SpecMTM(stats::ts(data))
  spec.mtm <- multitaper::spec.mtm(stats::ts(data), plot = FALSE)

  expect_true(utils::hasName(actual, "dof"))

  actual$dof <- NULL

  expect_equal(names(actual), names(spec.mtm))

})

test_that("log-smoothing works", {

  # only simple testing here on correct output structure

  data <- rnorm(1000)

  spec <- SpecMTM(stats::ts(data))
  actual <- LogSmooth(spec)

  expect_type(actual, "list")
  expect_s3_class(actual, "spec")

  expect_true(all(utils::hasName(actual, c("freq", "spec", "dof"))))

  expect_equal(spec$freq, actual$freq)

})

test_that("averaging spectra works", {

  m <- "MeanSpectrum: Spectra are of different lengths."
  s1 <- SpecMTM(stats::ts(rnorm(100)))
  s2 <- SpecMTM(stats::ts(rnorm(1000)))
  s3 <- SpecMTM(stats::ts(rnorm(1000)))

  expect_error(MeanSpectrum(list(s1, s2, s3)), m, fixed = TRUE)

  s1 <- SpecMTM(stats::ts(rnorm(1000)))
  s2 <- SpecMTM(stats::ts(rnorm(1000)))
  s3 <- SpecMTM(stats::ts(rnorm(1000)))

  m <- MeanSpectrum(list(s1, s2, s3))

  expect_type(m, "list")
  expect_s3_class(m, "spec")

  expect_true(all(utils::hasName(m, c("freq", "spec", "dof"))))

  expect_equal(m$freq, s1$freq)
  expect_equal(m$spec, (s1$spec + s2$spec + s3$spec) / 3)
  expect_equal(m$dof, s1$dof + s2$dof + s3$dof)

})

test_that("interpolating spectra works", {

  target <- list(freq = seq(0.1, 0.5, 0.1))

  x <- list(freq = c(0.1, 0.2, 0.4, 0.5, 0.6),
            spec = c(1, 2, 4, 5, 6))

  actual <- InterpolateSpectrum(x, target)

  expect_type(actual, "list")
  expect_s3_class(actual, "spec")
  expect_true(all(utils::hasName(actual, c("freq", "spec"))))

  expect_equal(actual$freq, target$freq)
  expect_equal(actual$spec, 1 : 5)

  x <- list(freq = c(0, 0.05, 0.1, 0.2, 0.4, 0.5),
            spec = c(pi, exp(1), 1, 2, 4, 5))

  actual <- InterpolateSpectrum(x, target)

  expect_equal(actual$freq, target$freq)
  expect_equal(actual$spec, 1 : 5)

  x <- list(freq = c(0.2, 0.4, 0.5, 0.6),
            spec = c(2, 4, 5, 6))

  expect_warning(
    actual <- InterpolateSpectrum(x, target),
    "NAs produced in interpolation as frequency axes do not overlap.",
    fixed = TRUE)

  expect_equal(actual$freq, target$freq)
  expect_equal(actual$spec, c(NA, 2 : 5))

})
