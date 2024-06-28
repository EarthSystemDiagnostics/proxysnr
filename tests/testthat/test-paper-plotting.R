test_that("paper plotting functions work", {

  data <- WrapSpectralResults(
    dml1 = dml$dml1, dml2 = dml$dml2, wais = wais,
    diffusion = diffusion.tf,
    time.uncertainty = time.uncertainty.tf,
    df.log = c(0.15, 0.15, 0.1)
  )

  # paper figure 2
  
  expect_no_error(
    # suppress warnings from neagtive values on log plot
    suppressWarnings(muench_laepple_fig02(data))
  )

  # paper figure 5

  expect_error(PublicationSNR(spec = data, data = "foobar"))

  corrupted_data <- data
  corrupted_data$dml1$corr.full <- NULL
  m <- "No version `corr.full` available for dataset `dml1`."
  expect_error(PublicationSNR(corrupted_data), m, fixed = TRUE)

  corrupted_data <- data
  corrupted_data$wais$raw <- NULL
  m <- "No version `raw` available for dataset `wais`."
  expect_error(PublicationSNR(corrupted_data, data = "raw"), m, fixed = TRUE)

  expect_no_error(
    # suppress warnings from permil sign conversion failures
    suppressWarnings(muench_laepple_fig05(proxysnr:::PublicationSNR(data)))
  )

})
