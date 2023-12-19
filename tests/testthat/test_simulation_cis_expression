# Test for correct output dimensions
test_that("output dimensions are correct", {
  n <- 1000
  p <- 40
  m <- 8
  G <- matrix(rbinom(n * p, size = 2, prob = 0.3), ncol = p)
  A <- matrix(sample(0:1, p * m, replace = TRUE), nrow = p)
  phi_v <- 0.05
  E <- simulation_cis_expression(G, A, phi_v)
  expect_equal(dim(E), c(n, m))
})

# Test for non-negative sigma2_cis
# Note: This test requires the function to return sigma2_cis values, which may require modification of the function
test_that("sigma2_cis values are non-negative", {
  n <- 1000
  p <- 40
  m <- 8
  G <- matrix(rbinom(n * p, size = 2, prob = 0.3), ncol = p)
  A <- matrix(sample(0:1, p * m, replace = TRUE), nrow = p)
  phi_v <- 0.05
  
  # Modify the function to return sigma2_cis values for testing
  E <- simulation_cis_expression(G, A, phi_v)
  # Assuming the modified function returns a list with E and sigma2_cis
  sigma2_cis_values <- E$sigma2_cis
  
  expect_true(all(sigma2_cis_values >= 0))
})

# Additional tests can be added as needed.

