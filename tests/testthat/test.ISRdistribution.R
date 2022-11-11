# Author: Quentin Grimonprez

context("ISR distribution")

test_that("simulISR samples mu when pi=1", {
  out <- simulISR(10, 1, 1:4)

  expect_equal(out, matrix(1:4, nrow = 10, ncol = 4, byrow = TRUE))
})

test_that("simulISR samples mu when pi=1", {
  set.seed(42)
  out <- simulISR(100000, 0.5, 1:3)
  x <- frequence(out)

  expect_lte(max(abs(x[, 4] / 100000 - rep(1 / 6, 6)) / rep(1 / 6, 6)), 2e-2)
})


test_that("probability works with a vector: univariate case", {
  # if pi == 1, the probability is 1 if x == mu, 0 otherwise
  p <- probability(1:4, 1:4, 1)
  expect_equal(p, 1)

  p <- probability(4:1, 1:4, 1)
  expect_equal(p, 0)

  p <- probability(1:4, 1:4, 0)
  expect_equal(p, 0)

  p <- probability(4:1, 1:4, 0)
  expect_equal(p, 1)
})

test_that("probability works with a vector: multivariate case", {
  # if pi == 1, the probability is 1 if x == mu, 0 otherwise
  p <- probability(c(1:4, 1:3), c(1:4, 1:3), c(1, 1), c(4, 3))
  expect_equal(p, 1)

  p <- probability(c(4:1, 1:3), c(1:4, 1:3), c(1, 1), c(4, 3))
  expect_equal(p, 0)
})

test_that("probability works with a matrix: univariate case", {
  # if pi == 1, the probability is 1 if x == mu, 0 otherwise
  p <- probability(matrix(1:4, nrow = 10, ncol = 4, byrow = TRUE), 1:4, 1)
  expect_equal(p, rep(1, 10))

  p <- probability(matrix(4:1, nrow = 10, ncol = 4, byrow = TRUE), 1:4, 1)
  expect_equal(p, rep(0, 10))

  # pi == 0.5 : uniform case
  x <- matrix(c(
    1, 2, 3,
    1, 3, 2,
    2, 1, 3,
    2, 3, 1,
    3, 1, 2,
    3, 2, 1
  ), ncol = 3, byrow = TRUE)

  mu <- 1:3
  p <- probability(x, mu, 0.5)
  expect_equal(p, rep(1 / nrow(x), nrow(x)))

  # pi = 0.9
  x <- matrix(c(
    1, 2, 3,
    1, 3, 2,
    2, 1, 3,
    2, 3, 1,
    3, 1, 2,
    3, 2, 1
  ), ncol = 3, byrow = TRUE)

  mu <- 1:3
  p <- probability(x, mu, 0.9)
  expect_equal(p, c(0.756, 0.084, 0.084, 0.036, 0.036, 0.004))
})


test_that("probability works with a matrix: multivariate case", {
  # if pi == 1, the probability is 1 if x == mu, 0 otherwise
  p <- probability(cbind(
    matrix(1:4, nrow = 10, ncol = 4, byrow = TRUE),
    matrix(1:3, nrow = 10, ncol = 3, byrow = TRUE)
  ), c(1:4, 1:3), c(1, 1), c(4, 3))
  expect_equal(p, rep(1, 10))

  p <- probability(cbind(
    matrix(4:1, nrow = 10, ncol = 4, byrow = TRUE),
    matrix(1:3, nrow = 10, ncol = 3, byrow = TRUE)
  ), c(1:4, 1:3), c(1, 1), c(4, 3))
  expect_equal(p, rep(0, 10))

  # pi == 0.5 : uniform case
  x <- matrix(c(
    1, 2, 3,
    1, 3, 2,
    2, 1, 3,
    2, 3, 1,
    3, 1, 2,
    3, 2, 1
  ), ncol = 3, byrow = TRUE)
  x <- cbind(x, x)

  mu <- 1:3
  p <- probability(x, c(mu, mu), c(0.5, 0.5), c(3, 3))
  expect_equal(p, rep(1 / (nrow(x) * nrow(x)), nrow(x)))
})
