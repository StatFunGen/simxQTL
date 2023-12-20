context("Simulate trans gene expression tests")

test_that("simulate_trans_expression handles T1E simulation correctly", {
    # Setup for T1E
    g <- 9
    n <- 1000
    A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
    phi_gene <- 0.15
    
    # Execute T1E simulation
    E <- simulate_trans_expression(A_trans, phi_gene, n = n, type = 'T1E')
    
    # Test if E has the correct dimensions
    expect_equal(dim(E), c(n, g))
})

test_that("simulate_trans_expression handles ACC/FDR simulation correctly", {
    # Setup for ACC/FDR
    m <- 8
    n <- 1000
    g <- 9
    E_cis <- matrix(rnorm(n * m), n, m)
    A_cis_trans <- matrix(sample(0:1, m * g, replace = TRUE), nrow = m)
    A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
    phi_gene <- 0.20
    
    # Execute ACC/FDR simulation
    E <- simulate_trans_expression(A_trans, phi_gene, E_cis
