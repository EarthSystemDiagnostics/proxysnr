test_that("CalculateTimeUncertaintyTF runs as expected", {

  skip_if_not_installed("simproxyage", minimum_version = "0.1.1")

  suppressMessages(
    actual <- CalculateTimeUncertaintyTF(nc = 3, acp = c(100, 50))
  )

  # OBS! To do the test here for nc = 3 is necessary since the function's
  # default of nc = 1 produces an error; see
  # <https://github.com/EarthSystemDiagnostics/simproxyage/issues/1>
  # Adjust this test once that big is fixed.
  # Also setting the acp manually here is necessary to avoid warnings, see same
  # issue.

  expect_type(actual, "list")
  expect_length(actual, 3)
  expect_named(actual, c("input", "stack", "ratio"))
  expect_true(is.spectrum(actual$input))
  expect_true(is.spectrum(actual$stack))
  expect_true(is.spectrum(actual$ratio))

  nt <- 100
  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))

  # test deprecated function name
  expect_warning(
    suppressMessages(
      actual.depr <- TimeUncertaintyTF(nc = 3, acp = c(100, 50))
    )
  )
  expect_equal(lapply(actual, "[[", "freq"), lapply(actual.depr, "[[", "freq"))

})
