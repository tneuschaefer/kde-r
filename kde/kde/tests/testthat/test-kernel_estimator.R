test_that("argument requirements work", {
  expect_error(kernel_estimator(x = TRUE))
  expect_error(kernel_estimator(x = numeric()))
  expect_error(kernel_estimator(x = NA, na.rm = FALSE))
  expect_error(kernel_estimator(x = 3, kernel = 4))
  expect_error(kernel_estimator(x = 3, kernel = stats::dnorm, bandwidth = TRUE))
  expect_error(kernel_estimator(x = 3, kernel = stats::dnorm, bandwidth = -1))
  expect_error(kernel_estimator(x = 3, kernel = stats::dnorm, bandwidth = numeric()))
})
