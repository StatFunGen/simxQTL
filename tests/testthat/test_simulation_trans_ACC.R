# Test for correct output dimensions
test_that("output dimensions are correct", {
  n <- 1000
  m <- 8
  g <- 9
  E_cis <- matrix(rnorm(n * m), n, m)
  A_cis_trans <- matrix(sample(0:1, m * g, replace = TRUE), nrow = m)
  A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
  phi_gene <- 0.20
  E_trans <- simulation_trans_expression_ACC_FDR(E_cis, A_cis_trans, A_trans, phi_gene)
  expect_equal(dim(E_trans), c(n, g))
})

# Additional tests can be added as needed.
