
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
simulate_cis_expression <- function(G, A, phi_v) {
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
        variance_sum <- var(G[, connected_snps, drop = FALSE] %*% beta[connected_snps])
        sigma2_cis <- var(G[, connected_snps[1]] * beta[connected_snps[1]]) / phi_v - variance_sum
        while (sigma2_cis <= 0) {
            phi_v <- phi_v - 0.01
            sigma2_cis <- var(E_comb[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum
        }
        
        # Simulate gene expression
        E_tmp <- G %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma2_cis))
        E[, j] <- scale(E_tmp)
    }
    
    return(E)
}



#' Simulate trans gene expression for different simulation scenarios
#'
#' This function can perform Type I error simulations or Accuracy and False Discovery Rate (ACC/FDR) simulations
#' depending on the specified type. For T1E simulations, it only needs an adjacency matrix and heritability. 
#' For ACC/FDR simulations, it requires additional cis gene expression and cis-trans adjacency matrix.
#'
#' @param A_trans Adjacency matrix among trans genes with dimensions g x g, required for both T1E and ACC/FDR.
#' @param phi_gene Per-gene heritability, required for both T1E and ACC/FDR.
#' @param n Sample size, required only for T1E simulations.
#' @param E_cis Matrix of cis gene expressions with dimensions n x m, required only for ACC/FDR simulations.
#' @param A_cis_trans Adjacency matrix between cis and trans genes with dimensions m x g, required for ACC/FDR simulations.
#' @param type Type of simulation to perform, either 'T1E' or 'ACC'.
#' @return A matrix E representing the simulated gene expression.
#'
#' @examples
#' # For T1E simulation
#' g <- 9
#' n <- 1000
#' A_trans <- matrix(sample(0:1, g * g, replace = TRUE), nrow = g)
#' phi_gene <- 0.15
#' E <- simulate_trans_expression(A_trans, phi_gene, n = n, type = 'T1E')
#'
#' # For ACC/FDR simulation
#' m <- 8
#' E_cis <- matrix(rnorm(n * m), n, m)
#' A_cis_trans <- matrix(sample(0:1, m * g, replace = TRUE), nrow = m)
#' E <- simulate_trans_expression(A_trans, phi_gene, E_cis = E_cis, A_cis_trans = A_cis_trans, type = 'ACC')
#' @export
simulate_trans_expression <- function(A_trans, phi_gene, n = NULL, E_cis = NULL, A_cis_trans = NULL, type = 'T1E') {
    if (type == 'T1E' && is.null(n)) {
        stop("For T1E simulation, 'n' must be provided.")
    }
    if (type == 'ACC' && (is.null(E_cis) || is.null(A_cis_trans))) {
        stop("For ACC/FDR simulation, 'E_cis' and 'A_cis_trans' must be provided.")
    }
    
    if (type == 'ACC') n <- nrow(E_cis)
    g <- ncol(A_trans)
    E <- if (type == 'T1E') matrix(rnorm(n * g), n, g) else matrix(0, n, g)
    
    if (type == 'T1E') {
        for (j in 1:g) {
            affecting_genes <- which(A_trans[, j] == 1)
            beta <- rep(0, g)
            if (length(affecting_genes) > 0) {
                beta[affecting_genes[1]] <- 1
                for (i in affecting_genes[-1]) {
                    beta[i] <- sqrt(var(E[, affecting_genes[1]]) / var(E[, i]))
                }
                variance_sum <- var(E[, affecting_genes, drop = FALSE] %*% beta[affecting_genes])
                sigma2_trans <- var(E[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum
                while (sigma2_trans <= 0) {
                    phi_gene <- phi_gene - 0.01
                    sigma2_trans <- var(E[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum
                }
                # Simulate gene expression
                E_tmp <- E %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma2_trans))
                E[, j] <- scale(E_tmp)
            }
        }
    } else if (type == 'ACC') {
        m <- ncol(E_cis)
        E_comb <- cbind(E_cis, E)
        A_comb <- rbind(A_cis_trans, A_trans)
        for (j in 1:g) {
            affecting_genes <- which(A_comb[, j] == 1)
            beta <- rep(0, m + g)
            if (length(affecting_genes) > 0) {
                beta[affecting_genes[1]] <- 1
                for (i in affecting_genes[-1]) {
                    beta[i] <- sqrt(var(E_comb[, affecting_genes[1]]) / var(E_comb[, i]))
                }
                variance_sum <- var(E_comb[, affecting_genes, drop = FALSE] %*% beta[affecting_genes])
                sigma2_gene <- var(E_comb[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum
                while (sigma2_gene <= 0) {
                    phi_gene <- phi_gene - 0.01
                    sigma2_gene <- var(E_comb[, affecting_genes[1]] * beta[affecting_genes[1]]) / phi_gene - variance_sum
                }
                # Multi-regression model for gene expression
                E_tmp <- E_comb %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma2_gene))
                E_comb[, m+j] <- scale(E_tmp)
            }
        }
        E <- E_comb[, (m + 1):(m + g)]
    }
    return(E)
}



#' Simulate trans gene expression for a mixture of cell types
#'
#' This function simulates trans gene expression considering different cell types,
#' each with its own gene-gene regulatory network.
#'
#' @param E_cis Matrix of cis gene expressions with dimensions n x m.
#' @param phi_gene Per-gene heritability.
#' @param A_cell_1 List containing adjacency matrices for cell-type 1 (A_cis_trans and A_trans).
#' @param A_cell_2 Optionally, a similar list for cell-type 2.
#' @param omega Optionally, a weight for combining cell types. If NULL, drawn from U(0.2, 0.8).
#' @return A list containing matrices E_trans, E_trans_cell_1, and E_trans_cell_2.
#' @examples
#' n <- 1000
#' m <- 8
#' g <- 9
#' phi_gene <- 0.20
#' E_cis <- matrix(rnorm(n * m), n, m)
#' A_cell_1 <- list(A_cis_trans = matrix(sample(0:1, m * g, replace = TRUE), nrow = m),
#'                 A_trans = matrix(sample(0:1, g * g, replace = TRUE), nrow = g))
#' result <- simulate_trans_mixture_celltype(E_cis, phi_gene, A_cell_1)
#' @export
simulate_trans_mixture_celltype <- function(E_cis, phi_gene, A_cell_1, 
                                            A_cell_2 = NULL, omega = NULL, 
                                            p = c(1/3, 2/3)) {
    
    if (is.null(A_cell_2)) {
        # Generate random A for cell type 2 if A_cell_2 is not provided
        A_cell_2 <- get_random_A(A_cell_1, p)
    }
    
    
    # Simulate for cell type 1 and 2
    E_trans_cell_1 <- simulate_trans_expression(A_trans = A_cell_1$A_trans, phi_gene = phi_gene,
                                                E_cis = E_cis,A_cis_trans = A_cell_1$A_cis_trans, type = 'ACC')
    E_trans_cell_2 <- simulate_trans_expression(A_trans = A_cell_2$A_trans, phi_gene = phi_gene,
                                                E_cis = E_cis,A_cis_trans = A_cell_2$A_cis_trans, type = 'ACC')
    
    # Determine omega if not provided
    if (is.null(omega)) {
        omega <- runif(1, 0.2, 0.8)
    }
    
    # Weighted sum
    E_trans <- omega * E_trans_cell_1 + (1 - omega) * E_trans_cell_2
    
    return(list(E_trans = E_trans, E_trans_cell_1 = E_trans_cell_1, E_trans_cell_2 = E_trans_cell_2))
}


#' Generate random adjacency matrices for trans gene associations
#'
#' This function creates random adjacency matrices for cis-trans and trans-trans effects
#' based on specified probabilities, as illustrated in Figure 7(B) and (C).
#'
#' @param A_cell_1 List containing fixed adjacency matrices for cell-type 1 (A_cis_trans and A_trans).
#' @param p Vector of probabilities for the associations.
#' @return A list containing random adjacency matrices for cell-type 2 (A_cis_trans and A_trans).
#'
#' @examples
#' A_cell_1 <- list(A_cis_trans = matrix(1, nrow = 8, ncol = 9),
#'                  A_trans = matrix(0, nrow = 9, ncol = 9))
#' p <- c(1/4, 1/2, 1/2) # Example probabilities
#' A_cell_2 <- get_random_A(A_cell_1, p)
#' @export
get_random_A <- function(A_cell_1, p = c(1/3, 2/3)) {
    m <- nrow(A_cell_1$A_cis_trans)
    g <- ncol(A_cell_1$A_cis_trans)
    
    # Assuming the first probability is for A_cis_trans and the second for A_trans
    A_cis_trans_random <- matrix(0, nrow = m, ncol = g)
    A_cis_trans_random[,1:4] <- matrix(rbinom(m * 4, size = 1, prob = p[1]), nrow = m)
    
    # For A_trans, considering a multi-layered approach as per Figure 7(C)
    A_trans_random <-matrix(0, nrow = g, ncol = g)
    A_trans_random[1:4,5:7] <- matrix(rbinom(4*3, size = 1, prob = p[2]), nrow = 4)
    A_trans_random[5:7,8:9] <- matrix(rbinom(3*2, size = 1, prob = p[2]), nrow = 3)
    
    # Ensure A_trans_random values are within 0 and 1 after layering probabilities
    A_trans_random[A_trans_random > 1] <- 1
    
    return(list(A_cis_trans = A_cis_trans_random, A_trans = A_trans_random))
}
