test_that("basic parallel processing works", {
  skip_on_cran()  # Skip on CRAN out of common courtesy

  # Setup test data
  set.seed(1234)
  x <- seq(-9, 9, length.out = 1000)
  y <- sin(x) + rnorm(1000, 0, 0.1)
  dat <- cbind(x, y)

  # Test parallel vs non-parallel results match
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))
  ## Ensure cluster is stopped even if test fails

  set.seed(1234)
  fit_parallel <- lgspline(cbind(dat[,'x']),
                          dat[,'y'],
                          cl = cl,
                          K = 2)

  set.seed(1234)
  fit_serial <- lgspline(cbind(dat[,'x']),
                        dat[,'y'],
                        K = 2)

  # Compare results
  expect_equal(fit_parallel$ytilde,
              fit_serial$ytilde,
              tolerance = 1e-5)

  # Test predictions match
  newx <- matrix(seq(-9, 9, length.out = 100))
  pred_parallel <- predict(fit_parallel, newx)
  pred_serial <- predict(fit_serial, newx)

  expect_equal(pred_parallel, pred_serial, tolerance = 1e-5)
})
