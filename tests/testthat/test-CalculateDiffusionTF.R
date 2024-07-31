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

  m <- "`df.log` must be of length 1 or `NULL`."
  expect_error(CalculateDiffusionTF(nt = nt, nc = nc, ns = 1, sigma = sigma,
                                    df.log = c(0.05, 0.05)), m, fixed = TRUE)

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

  # check attributes
  actualAttr <- attributes(actual$signal)
  expect_named(actualAttr,
               c("names", "class", "version", "N.sim", "log-smooth"))
  expect_true(class(actual$ratio) == "spec")
  expect_true(startsWith(actualAttr$version, "Creation date:"))
  expect_equal(actualAttr$N.sim, "Number of simulations used: N = 5.")
  expect_equal(actualAttr$`log-smooth`, "Log-smooth applied: No.")

  expect_equal(actualAttr, attributes(actual$diffused))
  expect_equal(actualAttr[c("class", "version", "N.sim", "log-smooth")],
               attributes(actual$ratio)[c("class", "version",
                                          "N.sim", "log-smooth")])

  # test deprecated function name
  expect_warning(
    actual.depr <- DiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma))
  expect_equal(lapply(actual, "[[", "freq"), lapply(actual.depr, "[[", "freq"))

  # test if coherent = TRUE
  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma,
                                 coherent = TRUE)

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))

  # test with smoothing
  actual <- CalculateDiffusionTF(nt = nt, nc = nc, ns = ns, sigma = sigma,
                                 df.log = 0.05)

  expect_equal(lengths(lapply(actual, "[[", "freq"), use.names = FALSE),
               rep(nt / 2, 3))
  expect_equal(attr(actual$ratio, "log-smooth"),
               "Log-smooth applied: Yes (df.log = 0.05).")

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
