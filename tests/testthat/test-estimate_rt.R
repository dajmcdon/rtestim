test_that("test when R crashes", {
  # Passed: k = 0 but x specified
  estimate_rt(1:10, degree = 0, lambda = 1, x = 1:10)

  ### Crashes: k = 0 & x unspecified
  # estimate_rt(1:10, degree = 0, lambda = 1)

  # Passed: k = 0 and whatever but x specified
  buildDx(10, 0, 1:10)
  buildDx(10, -1, 1:10)
  buildDx(10, 1, 1:10)
  buildDx(10, 2, 1:10)

  # Crashes: k = whatever, x unspecified
  # buildDx(10, -1, NULL)
  # buildDx(10, 0, NULL)
  # buildDx(10, 3, NULL)
  # buildDx(10, 3, c())

  # Passed: k = whatever but x = double(0)
  buildDx(10, 1, double(0))
  buildDx(10, 0, double(0))
  buildDx(10, -1, double(0))
})
