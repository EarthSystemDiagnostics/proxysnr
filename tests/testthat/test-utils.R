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

})
