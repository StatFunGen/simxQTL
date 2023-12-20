
context("Testing simulate_trans_mixture_celltype function")

test_that("simulate_trans_mixture_celltype returns correctly structured output", {
  n <- 1000
  m <- 8
  g <- 9
  phi_gene <- 0.20
  E_cis <- matrix(rnorm(n * m), nrow = n, ncol = m)
  A_cell_1 <- list(
    A_cis_trans = matrix(sample(0:1, m * g, replace = TRUE), nrow = m),
    A_trans = matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
  )
  result <- simulate_trans_mixture_celltype(E_cis, phi_gene, A_cell_1)
  
  expect_is(result, "list", info = "The output should be a list.")
  expect_true(all(sapply(result, function(x) dim(x) == c(n, g))), info = "Each matrix in the list should have dimensions n by g.")
})

test_that("simulate_trans_mixture_celltype handles omega correctly", {
  n <- 10
  m <- 5
  g <- 4
  phi_gene <- 0.20
  E_cis <- matrix(rnorm(n * m), nrow = n, ncol = m)
  A_cell_1 <- list(
    A_cis_trans = matrix(sample(0:1, m * g, replace = TRUE), nrow = m),
    A_trans = matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
  )
  omega <- 0.5
  result_with_omega <- simulate_trans_mixture_celltype(E_cis, phi_gene, A_cell_1, omega = omega)
  
  expect_equal(result_with_omega$E_trans, omega * result_with_omega$E_trans_cell_1 + (1 - omega) * result_with_omega$E_trans_cell_2, 
               info = "E_trans should be a weighted sum of E_trans_cell_1 and E_trans_cell_2 using omega.")
})

test_that("get_random_A returns matrices with correct structure and probabilities", {
  m <- 8
  g <- 9
  p <- c(1/4, 1/2) # Probabilities for A_cis_trans and A_trans
  
  A_cell_1 <- list(
    A_cis_trans = matrix(1, nrow = m, ncol = g),
    A_trans = matrix(0, nrow = g, ncol = g)
  )
  
  A_cell_2 <- get_random_A(A_cell_1, p)
  
  expect_equal(dim(A_cell_2$A_cis_trans), c(m, g), info = "Dimension of A_cis_trans does not match.")
  expect_equal(dim(A_cell_2$A_trans), c(g, g), info = "Dimension of A_trans does not match.")
  
  # Expect the proportion of 1s to be approximately equal to the set probability within a tolerance
  expect_equal(mean(A_cell_2$A_cis_trans), p[1], tolerance = 0.05, 
               info = "Proportion of 1s in A_cis_trans does not match the probability p[1].")
  expect_equal(mean(A_cell_2$A_trans), p[2], tolerance = 0.05, 
               info = "Proportion of 1s in A_trans does not match the probability p[2].")
})
