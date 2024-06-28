test_that("CalculateDiffusionTF error checks work", {

  nt <- 100
  nc <- 2

  m <- "Invalid length of supplied vector of diffusion lengths."
  sigma <- 1 : 10
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma),
               m, fixed = TRUE)

  m <- "Invalid dimensions of supplied array of diffusion lengths."
  sigma <- matrix(1 : 12, 3, 4)
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma),
               m, fixed = TRUE)

  sigma <- seq(0, 8, length.out = 100)

  m <- "`window` must have length two."
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    window = 11), m, fixed = TRUE)
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    window = c(11, 13, 42)), m, fixed = TRUE)

  m <- "`window[2]` must be > `window[1]`."
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    window = c(11, 3)), m, fixed = TRUE)

  m <- "`window[1]` must be >= 1."
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    window = c(0, 3)), m, fixed = TRUE)

  m <- "`window[2]` is > total number of time points."
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    window = c(11, 145)), m, fixed = TRUE)

})

test_that("CalculateDiffusionTF runs as expected", {

  sigma <- proxysnr:::diffusion.length$dml1$B32

  nt <- length(sigma)
  nc <- 2
  ns <- 5

  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma)

  expect_type(actual, "list")
  expect_length(actual, 3)
  expect_named(actual, c("signal", "diffused", "ratio"))
  expect_true(is.spectrum(actual$signal))
  expect_true(is.spectrum(actual$diffused))
  expect_true(is.spectrum(actual$ratio))

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))

  # test deprecated function name
  expect_warning(
    actual.depr <- DiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma))
  expect_equal(lapply(actual, "[[", "freq"), lapply(actual.depr, "[[", "freq"))

  # test if coherent = TRUE
  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma,
                                 coherent = TRUE)

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))

  # test with subsetting window set
  window <- c(11, 110)
  n <- diff(window) + 1
  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma,
                                 window = window)

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(n / 2, 3))

  # test for only one core, nc = 1
  nc <- 1
  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma)

  expect_type(actual, "list")
  expect_length(actual, 3)
  expect_named(actual, c("signal", "diffused", "ratio"))
  expect_true(is.spectrum(actual$signal))
  expect_true(is.spectrum(actual$diffused))
  expect_true(is.spectrum(actual$ratio))

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))

})
