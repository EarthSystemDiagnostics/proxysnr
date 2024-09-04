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
