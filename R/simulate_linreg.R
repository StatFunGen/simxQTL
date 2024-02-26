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
#' parse_num_causal_snps(10)       # 10 SNPs
#' parse_num_causal_snps(“10”)     # can also acept string format
#' parse_num_causal_snps("50pct")    # 50% of observed SNPs
#' parse_num_causal_snps("5avg")     # Average of 5 SNPs, sampled from truncated Poisson distribution
#' parse_num_causal_snps("0avg")     # Invalid: Average number of causal SNPs must be at least 1
#' parse_num_causal_snps("101pct")   # Invalid: Percentage must be in (0, 100]
#' parse_num_causal_snps("invalid")  # Invalid format
#' @export 
parse_num_causal_snps <- function(value) {
  is_pct <- FALSE
  if (grepl("[^0-9]+(pct|avg)?$", value, ignore.case = TRUE)) {
    num_tmp <- as.numeric(gsub("[^0-9]+", "", value))
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
    num_tmp = as.numeric(gsub("[^0-9]+", "", value))
  }

  return(list(value = num_tmp, is_pct = is_pct))
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
#' @param G Genotype matrix
#' @param ncausal Output from function parse_num_causal_snps, how many variants have non-negative effects (being causal)
#' @param ntrait Number of simulated phenotypes (traits)
#' @param shared_pattern: if is "all", all traits will have the same causal variant(s) with non-zero effect. if is "random", all traits will have independent (random) causal variant(s)
#' @return Vector of causal effects.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' # Create a genotype matrix
#' G = matrix(rbinom(1000, 2, 0.5), nrow = 1000, ncol = 50) 
#' B = sim_beta(G, ncausal = 5, ntrait = 3, is_h2g_total = T, shared_pattern = "all")
#' B = sim_beta(G, ncausal = 1, ntrait = 5, is_h2g_total = F, shared_pattern = "random")
sim_beta = function(G, ncausal, ntrait = 1, is_h2g_total = TRUE, shared_pattern = "all"){
  n_snps = ncol(G)
  ncausal = parse_num_causal_snps(ncausal)
  n_causal = if(ncausal$is_pct){
    max(1, round(ncausal$value * n_snps))
  }else{
    min(ncausal$value, n_snps)
  }
  B = matrix(0, nrow =  ncol(G), ncol = ntrait)
  
  if(shared_pattern == "all"){
    if(is_h2g_total){
      causal_index = sample(seq_len(n_snps), n_causal)
      beta = numeric(n_snps)
      beta[causal_index] = 1
      for(i in 1:ntrait){
        B[,i] = beta
      }
    }else{
      causal_index = sample(seq_len(n_snps), n_causal)
    beta = numeric(n_snps)
    beta[causal_index[1]] = 1
    var_vector = apply(as.matrix(G[,causal_index]), 2, var)
    beta[causal_index] = sqrt(beta[causal_index[1]]^2 * var_vector[1] / var_vector)
    for(i in 1:ntrait){
        B[,i] = beta
      }
    }
    
  }else if(shared_pattern == "random"){
    if(is_h2g_total){
      for(i in 1:ntrait){
        causal_index = sample(seq_len(n_snps), n_causal)
        beta = numeric(n_snps)
        beta[causal_index] = 1

        B[,i] = beta
      }
    }else{
          for(i in 1:ntrait){
        causal_index = sample(seq_len(n_snps), n_causal)
      beta = numeric(n_snps)
      beta[causal_index[1]] = 1
      var_vector = apply(as.matrix(G[,causal_index]), 2, var)
      beta[causal_index] = sqrt(beta[causal_index[1]]^2 * var_vector[1] / var_vector)

        B[,i] = beta
      }
    }    
    
    
  }else{
    stop('Shared pattern must be "all" or "random"!')
  }
  return(B)
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
#' @param residual_corr Matrix of residual correlations (NULL for independent samples).
#' @return A list containing the simulated phenotypes matrix (t * n, t = trait number, n = sample size) (`P`) and residual variance (`residual_var`).
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' G = matrix(rbinom(1000, 2, 0.5), nrow = 1000, ncol = 50) 
#' # Simulating effect sizes for two traits
#' B = sim_beta(G, ncausal = 5, ntrait = 3, is_h2g_total = F, shared_pattern = "all")
#' P = sim_multi_traits(G, B, h2g = 0.1, is_h2g_total = T, max_h2g = 1)
sim_multi_traits = function(G, B, h2g, is_h2g_total = TRUE, max_h2g = 1,  residual_corr = NULL){
  if (!is_h2g_total) {
        max_causal <- max(apply(B, 2, function(x) sum(x != 0)))
        h2g <- min(h2g, max_h2g/max_causal)
    }                             
  P = matrix(0, nrow = ncol(B), ncol = nrow(G)) # row: traits, column: different subjects
  mu = G %*% B 
  sigma = numeric(length = ncol(B))
  for(i in 1:ncol(mu)){
  if(is_h2g_total){
    sigma[i] = sqrt(var(mu[,i]) * (1-h2g) / h2g)
  }else{
    first_index = min(which(B[,i]==1))
    if(var(G[,first_index])/h2g - var(mu[,i]) >=0){
    sigma[i] =  sqrt(var(G[,first_index])/h2g - var(mu[,i]))
    }else{
      stop("Per SNP heritability too large, residual variance will be less than 0.")
    }
    }
  }
  if(is.null(residual_corr)){
    residual_corr <- diag(length(sigma))
  } 
  residual_var <- sweep(sweep(residual_corr, 2, sigma, "*"), 1, sigma, "*")
  P = mu + mvrnorm(n = nrow(G), mu = rep(0, ncol(residual_var)), Sigma = residual_var)
  return(list(P = t(P), residual_var = residual_var))
}