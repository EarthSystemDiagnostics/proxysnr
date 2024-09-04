test_that("plotting array spectra works", {

  # test error checking

  m <- "`spec` must be a list."
  expect_error(PlotArraySpectra(spec = 1), m, fixed = TRUE)

  m <- "`spec` must have elements `single`, `mean` and `stack`."
  expect_error(PlotArraySpectra(list(foo = 1, bar = 1)), m, fixed = TRUE)
  expect_error(PlotArraySpectra(list(single = 1, mean = 1 : 10)),
               m, fixed = TRUE)

  mean <- stack <- list(freq = 1 : 10, spec = 1 : 10)
  m <- "`spec$mean` must be a list with elements `freq` and `spec` of equal length."
  expect_error(PlotArraySpectra(list(single = 1, mean = 1, stack = stack)),
               m, fixed = TRUE)
  m <- "`spec$stack` must be a list with elements `freq` and `spec` of equal length."
  expect_error(PlotArraySpectra(list(single = 1, mean = mean, stack = 1)),
               m, fixed = TRUE)

  m <- "`spec$single` must be a list of spectra."
  expect_error(PlotArraySpectra(list(single = 1, mean = mean, stack = stack)),
               m, fixed = TRUE)

  single <- list(a = list(freq = 1 : 10, spec = 1 : 10), b = "ufp",
                 c = list(freq = 1 : 10, spec = 1 : 10))
  m <- "Cannot plot `spec$single[[2]]`: no spectral object."
  expect_error(
    PlotArraySpectra(list(single = single, mean = mean, stack = stack)),
    m, fixed = TRUE)

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

})

test_that("plotting SNR spectra works", {

  # test error checking

  m <- "`spec` must be a list."
  expect_error(PlotSNR(spec = 1), m, fixed = TRUE)

  m <- "`spec[[1]]` must be a list."
  spec <- list(1, list(2))
  expect_error(PlotSNR(spec), m, fixed = TRUE)

  m <- "`spec[[1]]` must have an element `snr`."
  spec <- list(list(foo = 1), list(snr = 2), list(spec = 3, snr = 4))
  expect_error(PlotSNR(spec), m, fixed = TRUE)

  m <- "Cannot plot `spec[[2]]$snr`: no spectral object."
  spec <- list(a = list(snr = list(freq = 1 : 3, spec = 1 : 3)),
               b = list(snr = list(5)))
  expect_error(PlotSNR(spec), m, fixed = TRUE)

  # create toy data
  n <- 100
  spec <- list(
    data1 = list(snr = list(freq = seq(0.01, 0.5, length.out = n),
                            spec = seq(1, 0.1, length.out = n))),
    data2 = list(snr = list(freq = seq(0.005, 0.5, length.out = n),
                            spec = seq(5, 0.1, length.out = n)))
  )

  expect_no_error(PlotSNR(spec))
  expect_no_error(PlotSNR(spec, col = 1))
  spec.no.names <- spec
  names(spec.no.names) <- NULL
  expect_no_error(PlotSNR(spec.no.names))

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

test_that("plotting stack correlation works", {

  # error checking

  msg <- "'data' must be a list."
  expect_error(PlotStackCorrelation(1), msg, fixed = TRUE)
  expect_error(PlotStackCorrelation(c(freq = "foo", correlation = "bar")),
               msg, fixed = TRUE)

  msg <- "'data' must have elements 'freq' and 'correlation'."
  expect_error(PlotStackCorrelation(list(foo = 1, bar = 2)), msg, fixed = TRUE)
  expect_error(PlotStackCorrelation(list(freq = 1, bar = 2)), msg, fixed = TRUE)

  freq <- seq(0.01, 0.5, length.out = 20)
  correlation <- matrix(runif(100), nrow = 5, ncol = 20)

  msg <- "Element 'correlation' must be a n x m matrix."
  expect_error(PlotStackCorrelation(list(freq = freq, correlation = 1 : 20)),
               msg, fixed = TRUE)

  msg <- paste("Length of element 'freq' must match",
               "number of columns of element 'correlation'.")
  expect_error(PlotStackCorrelation(
    list(freq = freq, correlation = matrix(1 : 20, nrow = 4, ncol = 5))),
    msg, fixed = TRUE)
  expect_error(PlotStackCorrelation(
    list(freq = 1 : 5, correlation = correlation)), msg, fixed = TRUE)

  data <- list(freq = freq, correlation = correlation)

  msg <- "Invalid x limit setting."
  expect_error(PlotStackCorrelation(data, xlim = 5), msg, fixed = TRUE)
  msg <- "Invalid y limit setting."
  expect_error(PlotStackCorrelation(data, ylim = 0.1), msg, fixed = TRUE)

  # test running function

  expect_no_error(PlotStackCorrelation(data))
  expect_no_error(PlotStackCorrelation(data, xlim = c(2, 5), ylim = c(10, 50)))
  expect_no_error(PlotStackCorrelation(data, xlim = c(1, NA), ylim = c(NA, 10)))

})

test_that("plotting tranfer functions works", {

  # test running function with defaults
  expect_no_error(PlotTF())

  # plot different number of datasets
  expect_no_error(PlotTF(dtf = diffusion.tf["dml1"],
                         ttf = time.uncertainty.tf[c("dml1", "wais")]))

  # plot with diffusion correction threshold
  expect_no_error(PlotTF(dtf.threshold = 0.67))

  # names must be created
  tf1 <- diffusion.tf
  tf2 <- time.uncertainty.tf
  names(tf1) <- NULL
  names(tf2) <- NULL
  expect_no_error(PlotTF(tf1, tf2))

  # names supplied as list
  expect_no_error(
    PlotTF(dtf = diffusion.tf["dml1"],
           ttf = time.uncertainty.tf[c("dml2", "wais")],
           names = list("bla", c("bloo", "blawp")))
  )

  # same number of datasets but names supplied as list -> still plot two legends
  expect_no_error(
    PlotTF(dtf = diffusion.tf[c("dml1", "wais")],
           ttf = time.uncertainty.tf[c("dml2", "wais")],
           names = list(c("Site A", "Site B"), c("Site C", "Site B")))
  )

  # name number mismatch
  msg1 <- "dtf: Number of data sets does not match number of names."
  expect_warning(PlotTF(dtf = diffusion.tf["dml1"], names = c("dml1", "wais")),
                 msg1, fixed = TRUE)
  msg2 <- "ttf: Number of data sets does not match number of names."
  expect_warning(PlotTF(ttf = time.uncertainty.tf,
                        names = c("spock", "sybok")),
                 msg2, fixed = TRUE)
  expect_warning(
    expect_warning(
      PlotTF(dtf = diffusion.tf["dml1"], ttf = time.uncertainty.tf,
             names = list(c("spock", "sybok"), c("dml1"))),
      msg1, fixed = TRUE),
    msg2, fixed = TRUE)

})
