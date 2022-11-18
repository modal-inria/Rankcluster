# Author: Quentin Grimonprez

context("Rank Manipulation")

test_that("convertSingleRank converts rank from ordering to ranking and vice-versa", {
  out <- convertSingleRank(1:4)
  expect_equal(out, 1:4)

  out <- convertSingleRank(c(2, 1, 4, 3))
  expect_equal(out, c(2, 1, 4, 3))

  out <- convertSingleRank(c(3, 1, 4, 2))
  expect_equal(out, c(2, 4, 1, 3))

  out <- convertSingleRank(c(2, 4, 1, 3))
  expect_equal(out, c(3, 1, 4, 2))
})

test_that("convertRank converts rank in vector", {
  out <- convertRank(1:4)
  expect_equal(out, 1:4)

  out <- convertRank(c(2, 1, 4, 3))
  expect_equal(out, c(2, 1, 4, 3))

  out <- convertRank(c(3, 1, 4, 2))
  expect_equal(out, c(2, 4, 1, 3))

  out <- convertRank(c(2, 4, 1, 3))
  expect_equal(out, c(3, 1, 4, 2))
})


test_that("convertRank converts rank in matrix", {
  x <- matrix(c(
    1, 2, 3, 4,
    2, 1, 4, 3,
    3, 1, 4, 2,
    2, 4, 1, 3
  ), ncol = 4, byrow = TRUE)

  out <- convertRank(x)
  expectedOut <- matrix(c(
    1:4, c(2, 1, 4, 3),
    c(2, 4, 1, 3), c(3, 1, 4, 2)
  ), ncol = 4, byrow = TRUE)
  expect_equal(out, expectedOut)
})


test_that("checkRank detects if the vector is a rank", {
  expect_true(checkRank(1:4))
  expect_false(checkRank(1:4, 5))
  expect_true(checkRank(c(2, 4, 1, 3)))
  expect_false(checkRank(c(4, 0, 1, 3)))
  expect_false(checkRank(c(5, 2, 1, 3)))
  expect_false(checkRank(c(4, NA, 1, 3)))
})


test_that("checkPartialRank detects if the vector with 0 or NA is a rank", {
  expect_true(checkPartialRank(1:4))
  expect_false(checkPartialRank(1:4, 5))
  expect_true(checkPartialRank(c(2, 4, 1, 3)))
  expect_true(checkPartialRank(c(4, 0, 1, 3)))
  expect_false(checkPartialRank(c(5, 2, 1, 3)))
  expect_false(checkPartialRank(c(4, 1, 1, 3)))
  expect_true(checkPartialRank(c(4, 1, NA, 3)))
})


test_that("checkTiePartialRank detects if the vector with 0, NA or ties is a rank", {
  expect_true(checkTiePartialRank(1:4))
  expect_false(checkTiePartialRank(1:4, 5))
  expect_true(checkTiePartialRank(c(2, 4, 1, 3)))
  expect_true(checkTiePartialRank(c(4, 0, 1, 3)))
  expect_true(checkTiePartialRank(c(4, NA, 1, 3)))
  expect_false(checkTiePartialRank(c(5, 2, 1, 3)))
  expect_true(checkTiePartialRank(c(4, 1, 1, 3)))
})


test_that("completeRank compeltes the rank when only 1 element is missing", {
  expect_equal(completeRank(1:4), 1:4)
  expect_equal(completeRank(c(1, 2, 3, 4, 0)), 1:5)
  expect_equal(completeRank(c(3, 5, 0, 1, 2)), c(3, 5, 4, 1, 2))
})


test_that("frequence works", {
  X <- matrix(c(
    rep(1:4, 5),
    rep(4:1, 3)
  ), ncol = 4, byrow = TRUE)
  out <- frequence(X)

  expectedOut <- matrix(c(
    1:4, 5,
    4:1, 3
  ), ncol = 5, byrow = TRUE)
  expect_equal(out, expectedOut)
})


test_that("unfrequence works", {
  X <- matrix(c(
    1:4, 5,
    4:1, 3
  ), ncol = 5, byrow = TRUE)
  out <- unfrequence(X)

  expectedOut <- matrix(c(
    rep(1:4, 5),
    rep(4:1, 3)
  ), ncol = 4, byrow = TRUE)

  expect_equal(out, expectedOut)
})
