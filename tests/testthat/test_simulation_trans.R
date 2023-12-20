context("Simulate trans gene expression tests")

test_that("simulate_trans_expression handles T1E simulation correctly", {
    g <- 9
    n <- 1000
    phi_gene <- 0.15
    A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
    
    E_t1e <- simulate_trans_expression(A_trans, phi_gene, n = n, type = 'T1E')
    
    expect_equal(dim(E_t1e), c(n, g), info = "T1E simulation should produce a matrix of dimensions n x g.")
    expect_true(all(colMeans(E_t1e) == 0, na.rm = TRUE), info = "T1E simulation should produce a scaled matrix with column means approximately 0.")
})

test_that("simulate_trans_expression handles ACC simulation correctly", {
    m <- 8
    n <- 1000
    g <- 9
    phi_gene <- 0.20
    E_cis <- matrix(rnorm(n * m), n, m)
    A_cis_trans <- matrix(sample(0:1, m * g, replace = TRUE), nrow = m)
    A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
    
    E_acc <- simulate_trans_expression(A_trans, phi_gene, E_cis = E_cis, A_cis_trans = A_cis_trans, type = 'ACC')
    
    expect_equal(dim(E_acc), c(n, g), info = "ACC simulation should produce a matrix of dimensions n x g.")
    expect_true(all(colMeans(E_acc) == 0, na.rm = TRUE), info = "ACC simulation should produce a scaled matrix with column means approximately 0.")
})
