test_that("obtaining correlation of stack works", {

  signal <- noise <- list(freq = 1 : 10, spec = 1 : 10)
  input <- list(signal = signal, noise = noise)

  N <- 1
  expected <- list(
    freq = 2 : 10,
    correlation = matrix(1 / sqrt(2), nrow = N, ncol = 9)
  )
  actual <- ObtainStackCorrelation(input)

  expect_equal(actual, expected)

  # test deprecated function name
  expect_warning(actual <- StackCorrelation(input))
  expect_equal(actual, expected)

  N <- c(1, 2, 4)
  expected <- list(
    freq = 2 : 10,
    correlation = matrix(c(rep(1 / sqrt(2), 9),
                           rep(1 / sqrt(3 / 2), 9),
                           rep(1 / sqrt(5 / 4), 9)),
                         byrow = TRUE, nrow = 3, ncol = 9)
  )
  actual <- ObtainStackCorrelation(input, N = N)

  expect_equal(actual, expected)


})
