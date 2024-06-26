test_that("general plotting functions work", {

  # plotting array spectra
  # ----------------------------------------------------------------------------

  # create artificial data
  nc <- 5
  nt <- 1000
  clim <- as.numeric(stats::arima.sim(model = list(ar = 0.7), n = nt))
  noise <- as.data.frame(replicate(nc, rnorm(n = nt)))

  proxy <- clim + noise

  spec <- ObtainArraySpectra(cores = proxy, df.log = 0.05)

  # check if plotting runs without error
  expect_no_error(PlotArraySpectra(spec))
  expect_no_error(
    PlotArraySpectra(spec, marker = 1 / c(50, 5),
                     xtm = c(50, 10, 5), xtl = c("a", "b", "c"),
                     ytm = c(0.01, 0.1, 1, 10),
                     ytl = c(0.01, 0.1, 1, 10))
  )

  # plotting SNR
  # ----------------------------------------------------------------------------

  # create toy data
  n <- 100
  spec <- list(
    data1 = list(snr = list(freq = seq(0.01, 0.5, length.out = n),
                            spec = seq(1, 0.1, length.out = n))),
    data2 = list(snr = list(freq = seq(0.005, 0.5, length.out = n),
                            spec = seq(5, 0.1, length.out = n)))
  )

  expect_no_error(PlotSNR(spec))

  m <- "Number of data sets does not match supplied number of names."
  expect_warning(PlotSNR(spec, names = "foo"), m, fixed = TRUE)

  # specify some cutoff frequencies
  i1 <- which.min(abs(spec$data1$snr$freq - 1 / 5))
  i2 <- which.min(abs(spec$data2$snr$freq - 1 / 20))

  # but both are missing
  m <- "`f.cut = TRUE` but cutoff frequency missing in input."
  expect_warning(PlotSNR(spec, f.cut = TRUE), m, fixed = TRUE)
  
  # only one is missing
  spec$data1$f.cutoff <- i1
  expect_warning(PlotSNR(spec, f.cut = TRUE), m, fixed = TRUE)

  # the other is missing
  spec$data1$f.cutoff <- NULL
  spec$data2$f.cutoff <- i2
  expect_warning(PlotSNR(spec, f.cut = TRUE), m, fixed = TRUE)

  # now correct
  spec$data1$f.cutoff <- i1
  spec$data2$f.cutoff <- i2
  expect_no_error(PlotSNR(spec, f.cut = TRUE))

  expect_no_error(
    PlotSNR(spec, xtm = c(500, 50, 5), xtl = c(500, 50, 5),
            ytm = c(0.05, 0.5, 5), ytl = c("very low", "low", "high"))
  )

})
