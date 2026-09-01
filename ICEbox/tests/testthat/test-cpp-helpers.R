test_that("matrix helpers return expected values", {
  x <- matrix(c(1, 2, 3, 4, 6, 8), nrow = 2, byrow = TRUE)

  expect_equal(rowMeans(rowCenter_cpp(x)), c(0, 0), tolerance = 1e-12)
  expect_equal(colSds_cpp(x), apply(x, 2, stats::sd), tolerance = 1e-12)
  expect_equal(melt_ice_curves_cpp(x), as.vector(t(x)))
})

test_that("numerical derivative handles a linear function", {
  grid <- seq(-2, 2, length.out = 9)
  curves <- rbind(3 * grid + 1, -2 * grid + 4)

  result <- derivative_cpp(curves, grid)

  expect_equal(result[1, ], rep(3, length(grid)), tolerance = 1e-10)
  expect_equal(result[2, ], rep(-2, length(grid)), tolerance = 1e-10)
})
