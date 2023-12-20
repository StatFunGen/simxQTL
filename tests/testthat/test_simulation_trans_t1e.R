# Test for correct output dimensions
test_that("output dimensions are correct", {
  g <- 9
  n <- 1000
  A <- matrix(0, nrow = g, ncol = g)
  A[5:9, 5:9] <- matrix(sample(0:1, 5 * 5, replace = TRUE), nrow = 5)
  phi_trans <- 0.15
  E <- simulation_trans_expression_t1e(A, phi_trans, n)
  expect_equal(dim(E), c(n, g))
})

# Additional tests can be added as needed.