
#' Simulate cis gene expression
#'
#' This function simulates cis gene expression based on genotype data.
#' It generates gene expression matrix E using genotype matrix G and adjacency matrix A.
#'
#' @param G A matrix of genotypes with dimensions n x p.
#' @param A A binary adjacency matrix of dimensions p x m indicating direct effects.
#' @param phi_v The per-SNP heritability value.
#'
#' @return A matrix E of dimensions n x m representing the simulated gene expression.
#'
#' @examples
#' n <- 1000
#' p <- 40
#' m <- 8
#' G <- matrix(rbinom(n * p, size = 2, prob = 0.3), ncol = p)
#' A <- matrix(sample(0:1, p * m, replace = TRUE), nrow = p)
#' phi_v <- 0.05
#' E <- simulation_cis_expression(G, A, phi_v)
#'
#' @export
simulation_cis_expression <- function(G, A, phi_v) {
  n <- nrow(G)
  m <- ncol(A)
  E <- matrix(NA, n, m)

  for (j in 1:m) {
    connected_snps <- which(A[, j] == 1)
    num_snps <- length(connected_snps)

    beta <- rep(0, ncol(G))
    for (i in connected_snps) {
      if (i == connected_snps[1]) {
        beta[i] <- 1
      } else {
        beta[i] <- sqrt(beta[connected_snps[1]]^2 * var(G[, connected_snps[1]]) / var(G[, i]))
      }
    }

    # Corrected calculation of variance_sum and sigma2_cis
    variance_sum <- var(G[, connected_snps] %*% beta[connected_snps])
    sigma2_cis <- var(G[, connected_snps[1]] * beta[connected_snps[1]]) / phi_v - variance_sum

    # Simulate gene expression
    E[, j] <- G %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma2_cis))
  }

  return(E)
}
