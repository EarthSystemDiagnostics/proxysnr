test_that("checking for spectral object works", {

  m <- "must be a list with elements `freq` and `spec` of equal length."

  expect_error(is.spectrum(1), paste("`1`", m), fixed = TRUE)
  expect_error(is.spectrum("char"), paste("`\"char\"`", m), fixed = TRUE)

  foo <- 1 : 10
  expect_error(is.spectrum(foo), paste("`foo`", m), fixed = TRUE)
  bar <- matrix(1 : 10, 5, 2)
  expect_error(is.spectrum(bar), paste("`bar`", m), fixed = TRUE)

  yet_another_bogus_var <- list(foo = 1 : 10)
  expect_error(is.spectrum(yet_another_bogus_var),
               paste("`yet_another_bogus_var`", m), fixed = TRUE)

  wrong_name <- list(freq = 1 : 10, specs = 1 : 10)
  expect_error(is.spectrum(wrong_name),
               paste("`wrong_name`", m), fixed = TRUE)

  different_lengths <- list(freq = 1 : 1000, spec = 1 : 10)
  expect_error(is.spectrum(different_lengths),
               paste("`different_lengths`", m), fixed = TRUE)

  s <- list(freq = 1 : 10, spec = 1 : 10)

  expect_true(is.spectrum(s))
  expect_true(is.spectrum(as.data.frame(s)))

})
