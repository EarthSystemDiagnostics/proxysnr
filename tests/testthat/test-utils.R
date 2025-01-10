test_that("checking if argument is a spectral object works", {

  m <- "must be a list with elements `freq` and `spec` of equal length."

  expect_error(check.if.spectrum(1), paste("`1`", m), fixed = TRUE)
  expect_error(check.if.spectrum("char"), paste("`\"char\"`", m), fixed = TRUE)

  foo <- 1 : 10
  expect_error(check.if.spectrum(foo), paste("`foo`", m), fixed = TRUE)
  bar <- matrix(1 : 10, 5, 2)
  expect_error(check.if.spectrum(bar), paste("`bar`", m), fixed = TRUE)

  yet_another_bogus_var <- list(foo = 1 : 10)
  expect_error(check.if.spectrum(yet_another_bogus_var),
               paste("`yet_another_bogus_var`", m), fixed = TRUE)

  wrong_name <- list(freq = 1 : 10, specs = 1 : 10)
  expect_error(check.if.spectrum(wrong_name),
               paste("`wrong_name`", m), fixed = TRUE)

  different_lengths <- list(freq = 1 : 1000, spec = 1 : 10)
  expect_error(check.if.spectrum(different_lengths),
               paste("`different_lengths`", m), fixed = TRUE)

  s <- list(freq = 1 : 10, spec = 1 : 10)

  expect_true(check.if.spectrum(s))
  expect_true(check.if.spectrum(as.data.frame(s)))

  # test on actual data

  spec <- ObtainArraySpectra(dml$dml1, df.log = 0.05)
  check.if.spectrum(spec$mean)
  check.if.spectrum(spec$stack)

  spec <- ObtainArraySpectra(dml$dml1)
  # <- without log-smooth, spec$stack is a more complex list with many elements
  # from spec.mtm function, so test needs to acount for this
  check.if.spectrum(spec$mean)
  check.if.spectrum(spec$stack)

})

test_that("testing for spectral object works", {

  expect_false(is.spectrum(1))
  expect_true(is.spectrum(list(freq = 1 : 10, spec = rnorm(10))))

})

test_that("checking of common frequency axes works", {

  target <- list(freq = seq(0.12, 0.2, 0.01))

  # target has lower bound than x
  expect_false(has.common.freq(x = list(freq = c(0.15, 0.3)), target))

  # target has higher bound than x
  expect_false(has.common.freq(x = list(freq = c(0.1, 0.15)), target))

  # target has both lower and higher bound than x
  expect_false(has.common.freq(x = list(freq = seq(0.15, 0.19, 0.001)), target))

  # target is completely outside range of x
  expect_false(has.common.freq(x = list(freq = seq(0.05, 0.1, 0.02)), target))
  expect_false(has.common.freq(x = list(freq = seq(0.3, 0.5, 0.1)), target))

  # matching axes
  expect_true(has.common.freq(x = list(freq = seq(0.12, 0.2, 0.005)), target))
  expect_true(has.common.freq(x = list(freq = seq(0, 0.5, 0.1)), target))
  expect_true(has.common.freq(x = list(freq = seq(0, 0.2, 0.001)), target))
  expect_true(has.common.freq(x = list(freq = seq(0.12, 1, 0.02)), target))

  # tests related to floating point representation issues
  spec <- ObtainArraySpectra(dml$dml2)
  tf <- time.uncertainty.tf[[2]]

  expect_true(has.common.freq(tf, spec$mean))

  spec <- ObtainArraySpectra(wais)
  tf <- time.uncertainty.tf[[3]]

  expect_true(has.common.freq(tf, spec$mean))
  # <- is FALSE for base >= type comparison due to floating point representation

})

test_that("subsetting spectrum by frequency range works", {

  s <- list(freq = seq(1 / 100, 0.5, length.out = 50), spec = rnorm(50))

  # per default no subset should be made
  expect_equal(fwindow(s), s)

  # use lower bound
  s1 <- list(freq = s$freq[-(1 : 20)], spec = s$spec[-(1 : 20)])
  expect_equal(fwindow(s, f.start = 0.21), s1)

  # use upper bound
  s2 <- list(freq = s$freq[1 : 42], spec = s$spec[1 : 42])
  expect_equal(fwindow(s, f.end = 0.42), s2)

  # use lower and upper bound
  s3 <- list(freq = s$freq[15 : 30], spec = s$spec[15 : 30])
  expect_equal(fwindow(s, f.start = 0.15, f.end = 0.3), s3)

})

test_that("fitting a power law to spectrum works", {

  x <- SimPLS(1000, beta = 1, alpha = 10) %>%
    SpecMTM(deltat = 1)

  plpar <- fit.powerlaw(x)

  expect_type(plpar, "list")
  expect_length(plpar, 2)
  expect_named(plpar, c("alpha", "beta"))
  expect_equal(lengths(plpar, use.names = FALSE), c(1, 1))

  expect_true(plpar$alpha > 0)
  expect_true(plpar$beta > 0)

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

  # test on case where floating point precision is an issue

  spec <- ObtainArraySpectra(wais)

  expected <- time.uncertainty.tf$wais
  actual   <- InterpolateSpectrum(time.uncertainty.tf$wais, spec$mean)

  expect_equal(actual$spec, expected$spec)

})

test_that("checking for array.par attribute works", {

  m <- "Attribute `array.par` missing from input object `spectra`."

  spectra  <- 1
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)
  attr(spectra, "array.foo") <- "bar"
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)

  m <- paste("Attribute `array.par` must a named vector",
             "with elements `nc`, `nt`, `res`.")

  attr(spectra, "array.par") <- "bar"
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)
  attr(spectra, "array.par") <- c(nt = 10)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)
  attr(spectra, "array.par") <- c(nc = 10, res = 1)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)

  m <- paste("Element `nc` of `array.par attribute",
             "(number of proxy records) must be a single integer >= 2.")
  attr(spectra, "array.par") <- c(nc = 0, nt = 100, res = 1)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)

  m <- paste("Element `nt` of `array.par attribute",
             "(number of observations per proxy record)",
             "must be a single integer > 8.")
  attr(spectra, "array.par") <- c(nc = 5, nt = 5, res = 1)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)
  attr(spectra, "array.par") <- c(nc = 5, nt = Inf, res = 1)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)

  m <- paste("Element `res` of `array.par attribute",
             "(resolution of proxy records)",
             "must be a single integer > 0.")
  attr(spectra, "array.par") <- c(nc = 5, nt = 100, res = 0)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)
  attr(spectra, "array.par") <- c(nc = 5, nt = 100, res = NA)
  expect_error(has.array.attribute(spectra), m, fixed = TRUE)

  attr(spectra, "array.par") <- c(nc = 5, nt = 100, res = 1)
  expect_no_error(has.array.attribute(spectra))
  expect_true(has.array.attribute(spectra))

})
