# The following is modified using chatgpt4 from
# https://github.com/mancusolab/twas_sim/blob/master/sim.py (Dec 16, 2023)

#' Parse Number of Causal SNPs
#'
#' Helper function to parse the number of causal SNPs with percentage and average modifiers.
#'
#' @param value Character string representing the number of causal SNPs.
#' @return A list containing the parsed number of causal SNPs and the type (percent or count).
#'
#' @examples
#' parse_num_causal_snps("10")       # 10 SNPs
#' parse_num_causal_snps("50pct")    # 50% of observed SNPs
#' parse_num_causal_snps("5avg")     # Average of 5 SNPs, sampled from truncated Poisson distribution
#' parse_num_causal_snps("0avg")     # Invalid: Average number of causal SNPs must be at least 1
#' parse_num_causal_snps("101pct")   # Invalid: Percentage must be in (0, 100]
#' parse_num_causal_snps("invalid")  # Invalid format
#' @export 
parse_num_causal_snps <- function(value) {
  is_pct <- FALSE
  if (grepl("^[0-9]+(pct|avg)?$", value, ignore.case = TRUE)) {
    num_tmp <- as.numeric(gsub("^[0-9]+", "", value))
    num_mod <- tolower(gsub("[^a-zA-Z]+", "", value))

    if (num_mod %in% c("pct", "avg")) {
      if (num_mod == "pct") {
        if (num_tmp <= 0 || num_tmp > 100) {
          stop("Percentage of causal SNPs must be in (0, 100].")
        }
        num_tmp <- num_tmp / 100
        is_pct <- TRUE
      } else if (num_mod == "avg") {
        if (num_tmp == 0) {
          stop("Average number of causal SNPs must be at least 1")
        }
        num_smpl <- 0
        while (num_smpl == 0) {
          num_smpl <- rpois(1, num_tmp)
        }
        num_tmp <- num_smpl
      }
    } else {
      if (num_tmp < 1) {
        stop("Number of causal SNPs must be at least 1")
      }
    }
  } else {
    stop("Invalid number of causal SNPs")
  }

  return(list(value = num_tmp, is_pct = is_pct))
}

#' Compute LD Matrix
#'
#' Computes the LD matrix from genotype data, adjusting for minor allele frequencies,
#' regularizing to ensure positive semi-definiteness, and standardizing to a correlation matrix.
#'
#' @param G Genotype matrix.
#' @param ld_ridge Ridge regularization parameter for LD matrix.
#' @return LD matrix.
#' @importFrom stats sd
#' @export
#'
#' @examples
#' # Random genotype matrix as an example
#' G <- matrix(rbinom(1000, 2, 0.5), nrow = 100, ncol = 10)
#' ld_ridge <- 0.1
#' LD <- get_ld(G, ld_ridge)
get_ld <- function(G, ld_ridge=0.1) {
  n <- nrow(G)
  p <- ncol(G)

  # Adjusting for minor allele frequencies
  mafs <- colMeans(G) / 2
  G <- sweep(G, 2, mafs * 2, FUN = "-")
  G <- sweep(G, 2, apply(G, 2, sd), FUN = "/")

  # Regularize LD matrix and adjust to correlation matrix
  LD <- crossprod(G) / n + diag(ld_ridge, p)
  LD <- LD / (1 + ld_ridge)

  return(LD)
}

#' Compute Lower Cholesky Decomposition
#'
#' Computes the lower Cholesky decomposition of the LD matrix.
#'
#' @param R The LD matrix for which the lower Cholesky decomposition is to be computed.
#' @return Lower Cholesky factor of the LD matrix.
#' @importFrom Matrix chol
#' @export
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)  # Example LD matrix
#' L <- get_lower_chol(R)                              # Compute lower Cholesky decomposition
get_lower_chol <- function(R) {
  L <- chol(R)
  return(L)
}

#' Compute Genetic Variance
#'
#' Compute genetic variance given betas and LD Cholesky factor.
#'
#' s2g := beta' V beta = beta' L L' beta
#'
#' @param RL Lower Cholesky factor of the p x p LD matrix for the population.
#' @param beta Genotype effects.
#' @return Genetic variance (s2g).
#' @export
#'
#' @examples
#' R <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)  # Example LD matrix
#' beta <- c(1, 2)                                 # Example genotype effects
#' RL <- get_lower_chol(R)                         # Compute lower Cholesky decomposition
#' compute_s2g(RL, beta)                            # Compute genetic variance
compute_s2g <- function(RL, beta) {
  RLtb <- crossprod(RL, beta) 
  s2g <- crossprod(RLtb)    
  return(s2g)
}

#' Simulate QTL Effects
#'
#' Sample QTL effects under a specified architecture.
#'
#' @param RL Lower Cholesky factor of the p x p LD matrix for the population.
#' @param ncausal Object of class NumCausalSNPs containing the number of causal SNPs to select.
#' @param eqtl_h2 The heritability of gene expression.
#' @param rescale Logical indicating whether to rescale effects such that var(b V b) = h2 (default is TRUE).
#' @return Vector of causal effects.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' # Create a correlation matrix as R
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
#' ncausal <- list(value = 50, is_pct = TRUE)  # 50% of observed SNPs
#' eqtl_h2 <- 0.5                              # Heritability of gene expression
#' RL <- get_lower_chol(R)                         # Compute lower Cholesky decomposition
#' sim_beta(RL, ncausal, eqtl_h2)               # Simulate QTL effects
sim_beta <- function(RL, ncausal, eqtl_h2, rescale = TRUE) {
  n_snps <- nrow(RL)
  
  n_qtls <- if (ncausal$is_pct) {
    max(1, round(ncausal$value * n_snps))
  } else {
    min(ncausal$value, n_snps)
  }
  
  if (eqtl_h2 != 0) {
    c_qtls <- sample(seq_len(n_snps), n_qtls)
    b_qtls <- numeric(n_snps)
    b_qtls[c_qtls] <- rnorm(n_qtls, mean = 0, sd = sqrt(eqtl_h2 / n_qtls))
    
    if (rescale) {
      s2g <- compute_s2g(RL, b_qtls)
      b_qtls <- b_qtls * sqrt(eqtl_h2 / s2g)
    }
  } else {
    b_qtls <- numeric(n_snps)
  }
  
  return(b_qtls)
}

#' Simulate Genotypes
#'
#' Sample genotypes from a multivariate normal approximation.
#'
#' @param RL Lower Cholesky factor of the p x p LD matrix for the population.
#' @param n Number of genotypes to sample.
#' @return Centered and scaled genotype matrix.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # LD matrix
#' n <- 100                                  # Number of genotypes to sample
#' RL <- get_lower_chol(R)                              # Compute lower Cholesky decomposition
#' G <- sim_geno(RL, n)                       # Simulate genotypes
sim_geno <- function(RL, n) {
  p <- nrow(RL)

  G <- RL %*% matrix(rnorm(n * p), ncol = p)
  G <- scale(G)  # Center and scale

  return(G)
}

#' Simulate a Polygenic Trait
#'
#' Simulate a complex trait as a function of latent genetic values and environmental noise.
#'
#' @param g Vector of latent genetic values.
#' @param h2g Heritability of the trait in the population.
#' @return Simulated phenotype.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # LD matrix
#' n <- 100                                  # Number of genotypes to sample
#' RL <- get_lower_chol(R)                              # Compute lower Cholesky decomposition
#' G <- sim_geno(RL, n)                       # Simulate genotypes
#' ncausal <- list(value = 50, is_pct = TRUE)  # 50% of observed SNPs
#' eqtl_h2 <- 0.5                              # Heritability of gene expression
#' b <- sim_beta(RL, ncausal, eqtl_h2) 
#' g <- G %*% b                   # latent genetic values
#' y <- simulate_polygenic_trait(g, eqtl_h2)  # Simulate the complex trait
simulate_polygenic_trait <- function(g, h2g) {
  n <- length(g)
  
  if (h2g > 0) {
    s2g <- var(g)  # Genetic variance
    s2e <- s2g * ((1 / h2g) - 1)  # Environmental variance
    e <- rnorm(n, mean = 0, sd = sqrt(s2e))
    y <- g + e
  } else {
    y <- rnorm(n, mean = 0, sd = 1)
  }
  
  # Standardize the phenotype
  y <- scale(y, center = TRUE, scale = FALSE)
  
  return(y)
}

#' Simulate GWAS Summary Statistics
#'
#' Simulate GWAS summary statistics directly using a multivariate normal approximation.
#' This method is efficient and designed for a large number of variants.
#'
#' @param RL Lower Cholesky factor of the LD matrix for the population.
#' @param ngwas Number of GWAS genotypes to sample.
#' @param beta Vector of latent eQTL effects for the causal gene.
#' @param h2ge Amount of phenotypic variance explained by the genetic component of gene expression.
#' @return A data frame containing estimated GWAS beta, standard error, and p-values.
#' @importFrom stats rnorm pnorm
#' @export
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Example LD matrix
#' ngwas <- 1000                             # Number of GWAS genotypes to sample
#' beta <- rnorm(2)                          # Latent eQTL effects
#' h2ge <- 0.5                               # Heritability of gene expression
#' RL <- get_lower_chol(R)                              # Compute lower Cholesky decomposition
#' sim_sumstats(RL, ngwas, beta, h2ge)        # Simulate GWAS summary statistics
sim_sumstats <- function(RL, ngwas, beta, h2ge) {
  n_snps <- nrow(RL)
  
  s2g <- sum((RL %*% beta)^2)
  
  if (h2ge > 0) {
    s2e <- s2g * ((1.0 / h2ge) - 1)
  } else {
    s2e <- 1.0
  }
  
  dof <- ngwas - 1
  tau2 <- s2e / ngwas
  se_gwas <- sqrt(rgamma(n_snps, shape = 0.5 * dof, rate = 0.5 * dof / tau2))
  DL <- se_gwas * RL
  
  beta_adj <- DL %*% t(RL) %*% solve(diag(se_gwas)) %*% beta
  b_gwas <- beta_adj + DL %*% matrix(rnorm(n_snps), ncol = 1)

  Z <- b_gwas / se_gwas
  pvals <- 2 * pnorm(abs(Z), lower.tail = FALSE)

  gwas <- data.frame(beta = b_gwas, se = se_gwas, pval = pvals)
  return(gwas)
}

# The following is modified using chatgpt4 from
# https://github.com/xueweic/fine-mcoloc/blob/b7aa9d27b0be13044224a4b2207717e74611d449/dsc/phenotype_simulation.r

#' Simulate Multiple Traits from Genotype and Effect Sizes
#'
#' This function simulates multiple traits (phenotypes) based on genotype data, 
#' a matrix of effect sizes, and heritability. It allows specifying if the heritability
#' is total or per-SNP, optionally scales the phenotypes, and can handle residual correlations.
#' per eQTL heritability is 0.05 to 0.07 according to https://www.nature.com/articles/ng.3506
#'
#' @param G Genotype matrix.
#' @param B Matrix of effect sizes for multiple traits.
#' @param h2g Heritability (proportion of variance explained by genetics).
#' @param is_h2g_total Logical indicating if h2g is total (TRUE) or per-SNP (FALSE).
#' @param max_h2g Maximum heritability allowed in the model (default is 0.8).
#' @param residual_corr Matrix of residual correlations (NULL for independent samples).
#' @param scale_Y Logical indicating if the phenotypes should be scaled (default is FALSE).
#' @return A list containing the simulated phenotypes (`Y`) and residual variance (`residual_var`).
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' G <- matrix(rnorm(1000), nrow = 100, ncol = 10) # Random genotype matrix
#' # Simulating effect sizes for two traits
#' B <- sapply(1:2, function(i) sim_beta(RL, ncausal, eqtl_h2)) 
#' h2g <- 0.5                                      # Heritability
#' result <- sim_multi_traits(G, B, h2g)           # Simulate multiple traits
sim_multi_traits <- function(G, B, h2g, is_h2g_total = TRUE, max_h2g = 0.8, 
                             residual_corr = NULL, scale_Y = FALSE) {
    if('scaled:scale' %in% names(attributes(G))) {
        SB <- sweep(B, 2, attr(G, "scaled:scale"), FUN = "/")
        GSB <- G %*% SB
        Yhat <- sweep(GSB, 2, colMeans(GSB), FUN = "-")
    } else {
        Yhat <- G %*% B
    }

    genetic_var <- apply(Yhat, 2, var)

    if (!is_h2g_total) {
        h2g <- min(h2g * ncol(B), max_h2g)
    }

    sigma <- sqrt(genetic_var / h2g - genetic_var)
    if (is.null(residual_corr)) {
        residual_corr <- diag(nrow(Yhat))
    }

    residual_var <- sweep(sweep(residual_corr, 2, sigma, "*"), 1, sigma, "*")
    Y <- Yhat + mvrnorm(n = nrow(G), mu = rep(0, ncol(residual_var)), Sigma = residual_var)

    if (scale_Y) {
        sd_Y <- apply(Y, 2, sd)
        Y <- sweep(Y, 2, sd_Y, FUN = "/")
        residual_var <- sweep(sweep(residual_var, 2, sd_Y, "/"), 1, sd_Y, "/")
        B <- sweep(B, 2, sd_Y, FUN = "/")
    }

    return(list(Y = Y, residual_var = residual_var))
}