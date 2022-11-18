# Author: Quentin Grimonprez

context("rankclust")

test_that("rankclust works on simulated data K=1", {
  set.seed(42)

  x <- simulISR(200, 0.8, 1:4)
  res <- rankclust(x)

  expect_s4_class(res, "Rankclust")
  expect_equal(res@K, 1)
  expect_equal(res@criterion, "bic")
  expect_true(res@convergence)
  expect_length(res@results, 1)
  expect_true(res@results[[1]]@convergence)
  expect_false(res@results[[1]]@partial)
  expect_equivalent(res@results[[1]]@pi, 0.8, tolerance = 1e-2)
  expect_equal(res@results[[1]]@proportion, 1)
  expect_equivalent(res@results[[1]]@mu, 1:4)
  expect_equal(res@results[[1]]@partition, rep(1, 200))
  expect_equal(res@results[[1]]@tik, matrix(1, nrow = 200))
})

test_that("rankclust works on simulated data K>1", {
  set.seed(42)

  x <- rbind(simulISR(200, 0.9, 1:4), simulISR(150, 0.9, c(2, 4, 3, 1)))
  res <- rankclust(x, K = 1:3)

  expect_s4_class(res, "Rankclust")
  expect_equal(res@K, 1:3)
  expect_equal(res@criterion, "bic")
  expect_true(res@convergence)
  expect_length(res@results, 3)
  expect_true(res@results[[2]]@convergence)
  expect_false(res@results[[2]]@partial)
  expect_equivalent(res@results[[2]]@mu, matrix(c(1:4, 2, 4, 3, 1), nrow = 2, byrow = TRUE))
})

test_that("rankclust finds the right parameters on a simple case", {
  set.seed(42)

  x <- rbind(simulISR(200, 0.95, 1:4), simulISR(150, 0.95, 4:1))
  res <- rankclust(x, K = 2)

  expect_s4_class(res, "Rankclust")
  expect_equal(res@K, 2)
  expect_equal(res@criterion, "bic")
  expect_true(res@convergence)
  expect_length(res@results, 1)
  expect_true(res@results[[1]]@convergence)
  expect_false(res@results[[1]]@partial)
  expect_equivalent(res@results[[1]]@mu, matrix(c(1:4, 4:1), nrow = 2, byrow = TRUE))
  expect_equal(res@results[[1]]@proportion, c(200, 150) / 350, tolerance = 2e-2)
  expect_equivalent(res@results[[1]]@pi, c(0.95, 0.95), tolerance = 2e-2)
  expect_gte(mean(res@results[[1]]@tik > 0.95) * 2, 0.9)
})


test_that("rankclust works on simulated data with tied", {
  set.seed(42)

  x <- simulISR(200, 0.9, 1:4)
  x[1, ] <- c(2, 1, 3, 3)
  x[2, ] <- c(1, 2, 2, 2)
  x[3, ] <- c(1, 1, 1, 1)
  x[4, ] <- c(1, 1, 3, 3)
  res <- rankclust(x, K = 1, criterion = "icl")

  expect_s4_class(res, "Rankclust")
  expect_equal(res@K, 1)
  expect_equal(res@criterion, "icl")
  expect_true(res@convergence)
  expect_length(res@results, 1)
  expect_true(res@results[[1]]@convergence)
  expect_true(res@results[[1]]@partial)
  expect_equivalent(res@results[[1]]@pi, 0.9, tolerance = 1e-2)
  expect_equal(res@results[[1]]@proportion, 1)
  expect_equivalent(res@results[[1]]@mu, 1:4)
  expect_equal(res@results[[1]]@partition, rep(1, 200))
  expect_equal(res@results[[1]]@tik, matrix(1, nrow = 200))
})


test_that("rankclust works on simulated data K>1", {
  set.seed(42)

  x <- rbind(simulISR(200, 0.9, 1:4), simulISR(150, 0.9, c(2, 4, 3, 1)))
  res <- rankclust(x, K = 1:3)

  expect_s4_class(res, "Rankclust")
  expect_equal(res@K, 1:3)
  expect_equal(res@criterion, "bic")
  expect_true(res@convergence)
  expect_length(res@results, 3)
  expect_true(res@results[[2]]@convergence)
  expect_false(res@results[[2]]@partial)
  expect_equivalent(res@results[[2]]@mu, matrix(c(1:4, 2, 4, 3, 1), nrow = 2, byrow = TRUE))
})
