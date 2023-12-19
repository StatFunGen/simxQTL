
# Testthat tests
# You can add these in a test file in the `tests/testthat` directory of your package

library(testthat)

test_that("simulate_causal_config returns correct output format", {
  A <- matrix(runif(100), ncol = 5)
  odds_ratio <- runif(5, 1, 2)
  baseline_odds <- 0.5
  n_samples <- 100
  result <- simulate_causal_config(A, odds_ratio, baseline_odds, n_samples)

  expect_true(is.list(result))
  expect_true(is.numeric(result$causal_probabilities))
  expect_true(is.matrix(result$samples))
  expect_equal(length(result$causal_probabilities), nrow(A))
  expect_equal(dim(result$samples), c(n_samples, nrow(A)))
})