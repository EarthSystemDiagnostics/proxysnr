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


})
