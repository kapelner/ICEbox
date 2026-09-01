test_that("ice builds curves and a PDP for a regression model", {
  X <- mtcars[c("wt", "hp", "disp")]
  fit <- stats::lm(mpg ~ wt + hp + disp, data = mtcars)

  result <- ice(
    fit,
    X = X,
    y = mtcars$mpg,
    predictor = "wt",
    num_grid_pts = 9,
    verbose = FALSE
  )

  expect_s3_class(result, "ice")
  expect_equal(dim(result$ice_curves), c(nrow(X), 9L))
  expect_equal(result$pdp, colMeans(result$ice_curves))
  expect_equal(result$xj, sort(X$wt))
  expect_equal(nrow(result$Xice), nrow(X))
})

test_that("ice validates predictors and incompatible transformations", {
  X <- mtcars[c("wt", "hp")]
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)

  expect_error(
    ice(fit, X, predictor = "missing", verbose = FALSE),
    "not found"
  )
  expect_error(
    ice(fit, X, predictor = "wt", logodds = TRUE, probit = TRUE, verbose = FALSE),
    "either logodds OR probit"
  )
})

test_that("dice estimates derivatives from ICE curves", {
  X <- mtcars[c("wt", "hp", "disp")]
  fit <- stats::lm(mpg ~ wt + hp + disp, data = mtcars)
  ice_result <- ice(
    fit,
    X = X,
    y = mtcars$mpg,
    predictor = "wt",
    verbose = FALSE
  )

  result <- dice(ice_result, sg_window_size = 5, verbose = FALSE)

  expect_s3_class(result, "dice")
  expect_equal(dim(result$d_ice_curves), c(nrow(X), length(ice_result$gridpts)))
  expect_length(result$dpdp, length(ice_result$gridpts))
  expect_length(result$actual_deriv, nrow(X))
})
