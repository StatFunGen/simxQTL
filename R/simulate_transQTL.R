
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


#' Simulate trans gene expression (Type I error simulations)
#'
#' This function simulates trans gene expression based on an adjacency matrix and per-gene heritability.
#' It generates an expression matrix E with standardized normal distribution for some genes,
#' and a multi-regression model for others.
#'
#' @param A The adjacency matrix of gene-gene regulatory causal relationship with dimensions g x g.
#' @param phi_gene The per-gene heritability value.
#' @param n The sample size.
#'
#' @return A matrix E of dimensions n x g representing the simulated gene expression.
#'
#' @examples
#' g <- 9
#' n <- 1000
#' A <- matrix(0, nrow = g, ncol = g)
#' A[5:9, 5:9] <- matrix(sample(0:1, 5 * 5, replace = TRUE), nrow = 5)
#' phi_gene <- 0.15
#' E <- simulation_trans_expression_t1e(A, phi_gene, n)
#'
#' @export
simulation_trans_expression_t1e <- function(A, phi_gene, n) {
  g <- ncol(A)
  E <- matrix(rnorm(n * g), n, g)  # Initialize E with standard normal distribution

  for (j in 1:g) {
    affecting_genes <- which(A[, j] == 1)
    beta <- rep(0, g)  # Initialize beta to zero

    if (length(affecting_genes) > 0) {
      beta[affecting_genes[1]] <- 1  # Set beta to 1 for the first affecting gene
      
      if (length(affecting_genes) > 1) {
        for (i in affecting_genes[-1]) {
          beta[i] <- sqrt(var(E[, affecting_genes[1]]) / var(E[, i]))
        }
      }

      # Calculate sigma^2_trans
      variance_sum <- var(E[, affecting_genes] %*% beta[affecting_genes])
      sigma2_trans <- var(E[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum

      # Multi-regression model for gene expression
      E[, j] <- rowSums(sapply(affecting_genes, function(i) E[, i] * beta[i])) + rnorm(n, mean = 0, sd = sqrt(sigma2_trans))
    }
  }

  return(E)
  

#' Simulate trans gene expression for accuracy and false discovery rate simulations
#'
#' This function simulates trans gene expression based on cis and trans gene expression data and adjacency matrices.
#' It creates a combined expression matrix and adjacency matrix for simulation.
#'
#' @param E_cis Matrix of cis gene expressions with dimensions n x m.
#' @param A_cis_trans Adjacency matrix between cis and trans genes with dimensions m x g.
#' @param A_trans Adjacency matrix among trans genes with dimensions g x g.
#' @param phi_gene Per-gene heritability.
#'
#' @return A matrix E_trans of dimensions n x g representing the simulated trans gene expression.
#'
#' @examples
#' n <- 1000
#' m <- 8
#' g <- 9
#' E_cis <- matrix(rnorm(n * m), n, m)
#' A_cis_trans <- matrix(sample(0:1, m * g, replace = TRUE), nrow = m)
#' A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
#' phi_gene <- 0.20
#' E_trans <- simulation_trans_expression_ACC_FDR(E_cis, A_cis_trans, A_trans, phi_gene)
#'
#' @export
simulation_trans_expression_ACC_FDR <- function(E_cis, A_cis_trans, A_trans, phi_gene) {
  n <- nrow(E_cis)
  m <- ncol(E_cis)
  g <- ncol(A_trans)
  E_trans <- matrix(0, n, g)  # Initialize E_trans with zeros
  E_comb <- cbind(E_cis, E_trans)  # Combine E_cis and E_trans

  # Create the block diagonal matrix A_comb
  A_comb_upper <- cbind(A_cis_trans, matrix(0, nrow = m, ncol = g))
  A_comb_lower <- cbind(matrix(0, nrow = g, ncol = m), A_trans)
  A_comb <- rbind(A_comb_upper, A_comb_lower)

  for (j in 1:(m + g)) {
    affecting_genes <- which(A_comb[, j] == 1)
    beta <- rep(0, m + g)  # Initialize beta to zero

    if (length(affecting_genes) > 0) {
      beta[affecting_genes[1]] <- 1  # Set beta to 1 for the first affecting gene
      
      if (length(affecting_genes) > 1) {
        for (i in affecting_genes[-1]) {
          beta[i] <- sqrt(var(E_comb[, affecting_genes[1]]) / var(E_comb[, i]))
        }
      }

      # Calculate sigma^2
      variance_sum <- var(E_comb[, affecting_genes] %*% beta[affecting_genes])
      sigma2_gene <- var(E_comb[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum

      # Multi-regression model for gene expression
      E_comb[, j] <- rowSums(sapply(affecting_genes, function(i) E_comb[, i] * beta[i])) + rnorm(n, mean = 0, sd = sqrt(sigma2_gene))
    }
  }

  # Extract the trans gene expression part from the combined matrix
  E_trans <- E_comb[, (m + 1):(m + g)]
  return(E_trans)
}
