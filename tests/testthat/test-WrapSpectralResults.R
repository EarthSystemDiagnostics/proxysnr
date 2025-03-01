test_that("WrapSpectralResults works", {

  # create artificial data
  
  nc <- c(3, 5, 11)
  nt <- c(200, 1000, 500)
  a1 <- c(0.7, 0.4, 0.85)

  makePseudoData <- function(i, nc, nt, a1) {

    clim <- as.numeric(stats::arima.sim(model = list(ar = a1[i]), n = nt[i]))
    noise <- as.data.frame(replicate(nc[i], rnorm(n = nt[i])))

    clim + noise

  }

  data <- lapply(seq_along(nc), makePseudoData, nc, nt, a1)

  # error checking

  m <- "No data sets supplied."
  expect_error(WrapSpectralResults(), m)

  m <- "Number of transfer functions must match number of datasets."
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   diffusion = NA),
               m, fixed = TRUE)
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   time.uncertainty = rep(NA, 4)),
               m, fixed = TRUE)
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   diffusion = list(1, 0.5),
                                   time.uncertainty = rep(NA, 4)),
               m, fixed = TRUE)

  m <- "Length of `df.log` must be 1 or equal to the number of datasets."
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   df.log = c(0.05, 0.1)), m, fixed = TRUE)
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   df.log = NULL), m, fixed = TRUE)
  m <- "Length of `res` must be 1 or equal to the number of datasets."
  expect_error(WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                                   res = 1 : 11), m, fixed = TRUE)

  # test without any correction

  spec <- WrapSpectralResults(data[[1]], data[[2]], data[[3]])

  expect_type(spec, "list")
  expect_length(spec, 3)
  expect_equal(length(names(spec)), 0)

  expect_equal(sapply(spec, names), rep("raw", 3))
  expect_equal(names(spec[[1]]$raw), c("signal", "noise", "snr", "f.cutoff"))
  expect_equal(names(spec[[2]]$raw), c("signal", "noise", "snr", "f.cutoff"))
  expect_equal(names(spec[[3]]$raw), c("signal", "noise", "snr", "f.cutoff"))

  expect_true(is.spectrum(spec[[1]]$raw$signal))
  expect_true(is.spectrum(spec[[1]]$raw$noise))
  expect_true(is.spectrum(spec[[1]]$raw$snr))
  expect_true(is.na(spec[[1]]$raw$f.cutoff))

  expect_equal(attr(spec[[1]]$raw, "array.par"),
               c(nc = nc[1], nt = nt[1], res = 1))
  expect_equal(attr(spec[[2]]$raw, "array.par"),
               c(nc = nc[2], nt = nt[2], res = 1))
  expect_equal(attr(spec[[3]]$raw, "array.par"),
               c(nc = nc[3], nt = nt[3], res = 1))

  # test with correction functions applied to one or to both datasets

  diffusion <- list(
    list(freq = seq(1 / nt[1], 0.5, length.out = nt[1] / 2),
         spec = rep(0.8, nt[1] / 2)),
    NA,
    list(freq = seq(1 / nt[3], 0.5, length.out = nt[3] / 2),
         spec = rep(0.8, nt[3] / 2))
  )
  time.uncertainty <- list(
    NA,
    list(freq = seq(1 / nt[2], 0.5, length.out = nt[2] / 2),
         spec = rep(0.95, nt[2] / 2)),
    list(freq = seq(1 / nt[3], 0.5, length.out = nt[3] / 2),
         spec = rep(0.95, nt[3] / 2))
  )

  spec <- WrapSpectralResults(data[[1]], data[[2]], data[[3]],
                              diffusion = diffusion,
                              time.uncertainty = time.uncertainty,
                              crit.diffusion = 1 / 0.9)

  expect_equal(names(spec[[1]]), c("raw", "corr.diff.only"))
  expect_equal(names(spec[[2]]), c("raw", "corr.t.unc.only"))
  expect_equal(names(spec[[3]]),
               c("raw", "corr.diff.only", "corr.t.unc.only", "corr.full"))

  expect_true(is.numeric(spec[[1]]$corr.diff.only$f.cutoff))
  expect_true(is.na(spec[[2]]$corr.t.unc.only$f.cutoff))
  expect_true(is.na(spec[[3]]$corr.t.unc.only$f.cutoff))
  expect_true(is.numeric(spec[[3]]$corr.diff.only$f.cutoff))
  expect_true(is.numeric(spec[[3]]$corr.full$f.cutoff))

})
