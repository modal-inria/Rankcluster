# Author: Quentin Grimonprez

context("Test functions")

test_that("kullback-leibler divergence is close to 0 for same parameters", {
  set.seed(42)
  proportion1 <- c(0.4, 0.6)
  pi1 <- matrix(c(0.8, 0.75), nrow = 2)
  mu1 <- matrix(c(1, 2, 3, 4, 4, 2, 1, 3), nrow = 2, byrow = TRUE)

  dK <- kullback(proportion1, pi1, mu1, proportion1, pi1, mu1, 4)
  expect_equal(dK, 0, tolerance = 1e-2)
})

test_that("kullback-leibler divergence is large for different parameters", {
  set.seed(42)
  proportion1 <- c(0.4, 0.6)
  pi1 <- matrix(c(0.8, 0.75), nrow = 2)
  mu1 <- matrix(c(
    1, 2, 3, 4,
    4, 2, 1, 3
  ), nrow = 2, byrow = TRUE)

  proportion2 <- c(0.8, 0.2)
  pi2 <- matrix(c(0.9, 0.9), nrow = 2)
  mu2 <- matrix(c(
    4, 3, 2, 1,
    1, 3, 4, 2
  ), nrow = 2, byrow = TRUE)

  dK <- kullback(proportion1, pi1, mu1, proportion2, pi2, mu2, 4)
  expect_gt(dK, 1.5)
})


test_that("kullback-leibler divergence is close to 0 for close parameters", {
  set.seed(42)
  proportion1 <- c(0.4, 0.6)
  pi1 <- matrix(c(0.8, 0.75), nrow = 2)
  mu1 <- matrix(c(
    1, 2, 3, 4,
    4, 2, 1, 3
  ), nrow = 2, byrow = TRUE)

  proportion2 <- c(0.42, 0.58)
  pi2 <- matrix(c(0.85, 0.7), nrow = 2)
  mu2 <- matrix(c(
    1, 2, 3, 4,
    4, 2, 1, 3
  ), nrow = 2, byrow = TRUE)

  dK <- kullback(proportion1, pi1, mu1, proportion2, pi2, mu2, 4)
  expect_lt(dK, 0.3)
})

test_that("khi2 test has a p-value greater than 0.05 for data from given parameters", {
  set.seed(42)
  proportion <- c(0.4, 0.6)
  pi <- c(0.8, 0.75)
  mu <- matrix(c(
    1, 2, 3, 4,
    4, 2, 1, 3
  ), nrow = 2, byrow = TRUE)

  data <- rbind(
    simulISR(proportion[1] * 500, pi[1], mu[1, ]),
    simulISR(proportion[2] * 500, pi[2], mu[2, ])
  )
  pval <- khi2(data, proportion, mu, pi)

  expect_gt(pval, 0.05)
})

test_that("khi2 test has a p-value lower than 0.05 for data from other parameters", {
  set.seed(42)
  proportion <- c(0.4, 0.6)
  pi <- c(0.8, 0.75)
  mu <- matrix(c(
    1, 2, 3, 4,
    4, 2, 1, 3
  ), nrow = 2, byrow = TRUE)

  data <- rbind(
    simulISR(0.8 * 500, 0.9, c(4, 3, 2, 1)),
    simulISR(0.2 * 500, 0.9, c(3, 1, 4, 2))
  )
  pval <- khi2(data, proportion, mu, pi)

  expect_lt(pval, 0.05)
})
