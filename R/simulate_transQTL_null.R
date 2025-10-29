#' Generate gene expression data for trans-QTL analysis
#'
#' This function simulates genotype and gene expression data for trans-QTL analysis, 
#' including both cis and trans effects with specified network structures.
#'
#' @param n Sample size
#' @param p Number of SNPs
#' @param A_snp_cis Binary adjacency matrix indicating SNP-cis gene relationships
#' @param A_cis_trans Binary adjacency matrix indicating cis-trans gene relationships
#' @param A_trans Binary adjacency matrix indicating trans-trans gene relationships
#' @param MAF Optional vector of minor allele frequencies. If NULL, MAFs will be generated randomly.
#' @param beta Effect size for SNP to cis gene effects
#' @param gamma_cis_trans Effect size for cis to trans gene effects
#' @param gamma_trans_trans Effect size for trans to trans gene effects
#' @param maf_setting Integer flag for MAF generation method (1=uniform random, 2=sample from provided MAF)
#'
#' @return A list containing:
#' \itemize{
#'   \item maf: Vector of minor allele frequencies
#'   \item G_unscaled: Unscaled genotype matrix
#'   \item G: Scaled genotype matrix
#'   \item G_centered: Column-centered genotype matrix
#'   \item E_cis_unscaled: Unscaled cis gene expression matrix
#'   \item E_cis: Scaled cis gene expression matrix
#'   \item E_cis_centered: Column-centered cis gene expression matrix
#'   \item E_trans_unscaled: Unscaled trans gene expression matrix
#'   \item E_trans: Scaled trans gene expression matrix
#'   \item E_trans_centered: Column-centered trans gene expression matrix
#' }
#'
#' @examples
#' # Define adjacency matrices
#' A_snp_cis <- as.matrix(Matrix::bdiag(lapply(1:8, function(x) rep(1,5))))
#' A_snp_cis <- rbind(A_snp_cis, matrix(0, nrow = 20, ncol = 8))
#' 
#' A_cis_trans <- cbind(c(1,1,0,0,0,0,0,0),
#'                     c(0,0,1,1,0,0,0,0),
#'                     c(0,0,0,0,1,1,0,0),
#'                     c(0,0,0,0,0,0,1,1),0,0,0,0,0)
#'
#' A_trans <- cbind(0,0,0,0,c(1,rep(0,8)),
#'                 c(0,1,1,rep(0,6)),
#'                 c(0,0,0,1,rep(0,5)),
#'                 c(0,0,0,0,1,1,0,0,0),
#'                 c(rep(0,5),1,1,0,0))
#'
#' # Simulate data
#' data <- gene_data(n = 1000, p = 60, A_snp_cis, A_cis_trans, A_trans,
#'                   beta = 0.5, gamma_cis_trans = 0.5, gamma_trans_trans = 0.5)
#'
#' @export
gene_data <- function(n, p, A_snp_cis, A_cis_trans, A_trans, MAF = NULL, 
                     beta = 1,
                     gamma_cis_trans = 0.5,
                     gamma_trans_trans = 0.5,
                     maf_setting = 1) {
    res.list <- list()
    
    # Generate SNPs based on MAF setting
    if(maf_setting == 1) {
        G <- lapply(1:p, function(x) rbinom(n, 2, runif(1, 0.1, 0.4)))
    } else {
        G <- lapply(1:p, function(x) {
            maf <- sample(MAF[MAF >= 0.01], 1)
            Gcol <- rbinom(n, 2, maf)
            while(sum(Gcol) == 0) {
                maf <- sample(MAF[MAF >= 0.01], 1)
                Gcol <- rbinom(n, 2, maf)
            }
            Gcol
        })
    }
    
    # Process genotype data
    G <- do.call(cbind, G)
    res.list$maf <- colMeans(G)/2
    res.list$G_unscaled <- G
    G_scaled <- scale(G)
    G_centered <- scale(G, scale = FALSE)
    res.list$G <- G_scaled
    res.list$G_centered <- G_centered
    
    # Generate cis effects
    E_cis <- matrix(nrow = n, ncol = ncol(A_snp_cis))
    for(j in 1:ncol(A_snp_cis)) {
        asso <- which(A_snp_cis[,j] == 1)
        beta_j <- rnorm(length(asso), mean = beta, sd = sqrt(0.04))
        mean_cis <- G[,asso, drop = FALSE] %*% as.matrix(beta_j)
        E_cis[,j] <- sapply(mean_cis, function(x) rnorm(1, mean = x, sd = 1))
    }
    res.list$E_cis_unscaled <- E_cis
    res.list$E_cis <- scale(E_cis)
    res.list$E_cis_centered <- scale(E_cis, scale = FALSE)
    
    # Generate trans effects - considering both cis and trans influences
    E_trans <- matrix(0, nrow = n, ncol = ncol(A_trans))
    for(j in 1:ncol(A_trans)) {
        asso_cis <- which(A_cis_trans[,j] == 1)
        asso_trans <- which(A_trans[,j] == 1)
        
        mean_total <- 0
        # cis to trans effects
        if(length(asso_cis) > 0) {
            beta_ct <- rnorm(length(asso_cis), mean = gamma_cis_trans, sd = 0)
            mean_total <- mean_total + E_cis[, asso_cis, drop = FALSE] %*% beta_ct
        }
        # trans to trans effects
        if(length(asso_trans) > 0) {
            beta_tt <- rnorm(length(asso_trans), mean = gamma_trans_trans, sd = 0)
            mean_total <- mean_total + E_trans[, asso_trans, drop = FALSE] %*% beta_tt
        }
        
        E_trans[,j] <- sapply(mean_total, function(x) rnorm(1, mean = x, sd = 1))
    }
    
    res.list$E_trans_unscaled <- E_trans
    res.list$E_trans <- scale(E_trans)
    res.list$E_trans_centered <- scale(E_trans, scale = FALSE)
    
    return(res.list)
}

#' Generate null data by random shuffling of genotypes
#'
#' This function creates a null dataset by randomly shuffling the genotype data
#' while maintaining the original gene expression data.
#'
#' @param data_h1 A dataset generated by the gene_data function
#'
#' @return A modified dataset with randomly shuffled genotypes
#'
#' @examples
#' # Generate original data
#' data_h1 <- gene_data(n = 1000, p = 60, A_snp_cis, A_cis_trans, A_trans)
#' 
#' # Generate null data by shuffling genotypes
#' data_null <- gene_data_null1(data_h1)
#'
#' @export
gene_data_null1 <- function(data_h1) {
    data <- data_h1
    data$G <- apply(data_h1$G, 2, sample)
    return(data)
}

#' Generate null data by removing cis-to-trans effects
#'
#' This function creates a null dataset by regenerating trans gene expression data
#' without the influence of cis genes, while maintaining trans-to-trans effects.
#'
#' @param data_h1 A dataset generated by the gene_data function
#' @param A_trans Binary adjacency matrix indicating trans-trans gene relationships
#' @param gamma_trans_trans Effect size for trans to trans gene effects
#'
#' @return A modified dataset with new trans gene expression lacking cis-to-trans effects
#'
#' @examples
#' # Generate original data
#' data_h1 <- gene_data(n = 1000, p = 60, A_snp_cis, A_cis_trans, A_trans)
#'
#' # Generate null data without cis-to-trans effects
#' data_null <- gene_data_null2(data_h1, A_trans, gamma_trans_trans = 0.5)
#'
#' @export
gene_data_null2 <- function(data_h1, A_trans, gamma_trans_trans = 0.5) {
    data <- data_h1
    n <- nrow(data_h1$E_trans)
    
    E_trans_null <- matrix(0, nrow = n, ncol = ncol(data_h1$E_trans))
    # Generate Genes 1-4 from standard normal
    E_trans_null[, 1:4] <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
    # Generate Genes 5-9 using only trans-gene relationships
    for(j in 5:9) {
        asso_trans <- which(A_trans[,j] == 1)
        if(length(asso_trans) > 0) {
            temp_trans <- E_trans_null[, asso_trans, drop = FALSE]
            beta <- rnorm(ncol(temp_trans), mean = gamma_trans_trans, sd = 0)
            mean_trans <- temp_trans %*% as.matrix(beta)
            E_trans_null[,j] <- sapply(mean_trans, function(x) rnorm(1, mean = x, sd = 1))
        } else {
            E_trans_null[,j] <- rnorm(n)
        }
    }
    data$E_trans <- scale(E_trans_null)
    return(data)
}

#' Generate null data using ARCHIE setting
#'
#' This function creates a null dataset by directly generating a Z-statistic matrix
#' under the multivariate normal assumption used in the ARCHIE method.
#'
#' @param n Sample size
#' @param p Number of SNPs
#' @param data_h1 A dataset generated by the gene_data function
#' @param A_snp_cis Binary adjacency matrix indicating SNP-cis gene relationships
#' @param A_cis_trans Binary adjacency matrix indicating cis-trans gene relationships
#' @param A_trans Binary adjacency matrix indicating trans-trans gene relationships
#'
#' @return A matrix of Z-statistics under the null hypothesis
#'
#' @examples
#' # Generate original data
#' data_h1 <- gene_data(n = 1000, p = 60, A_snp_cis, A_cis_trans, A_trans)
#' 
#' # Generate null Z-matrix using ARCHIE setting
#' Z_matrix_null <- gene_data_null3(1000, 60, data_h1, A_snp_cis, A_cis_trans, A_trans)
#'
#' @importFrom MASS mvrnorm
#' @export
gene_data_null3 <- function(n, p, data_h1, A_snp_cis, A_cis_trans, A_trans) {
    Sigma_EE <- cor(data_h1$E_trans)
    Sigma_GG <- diag(p)  
    
    Sigma_kronecker <- kronecker(Sigma_EE, Sigma_GG)
    vec_Sigma_GE <- MASS::mvrnorm(1, mu = rep(0, p * ncol(data_h1$E_trans)), 
                              Sigma = Sigma_kronecker)
    
    Z_matrix_null <- matrix(vec_Sigma_GE, nrow = p)
    return(Z_matrix_null)
}