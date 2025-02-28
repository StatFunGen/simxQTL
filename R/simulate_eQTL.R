###############################################################################
#' Simulate Sparse Effects for eQTL Data
#'
#' This function simulates the sparse component of the genetic effects,
#' including a sentinel SNP and a specified number of additional sparse SNPs.
#'
#' @param G A standardized genotype matrix.
#' @param h2_sparse Heritability allocated to sparse effects.
#' @param prop_h2_sentinel Proportion of h2_sparse attributed to the sentinel SNP.
#' @param n_other_sparse Number of additional sparse SNPs.
#' @return A list containing:
#'   \item{beta}{A vector of effect sizes with nonzero entries for the sparse SNPs.}
#'   \item{sentinel_index}{Index of the sentinel SNP.}
#'   \item{other_sparse_indices}{Indices of the additional sparse SNPs.}
#' @export
simulate_sparse_effects <- function(G, h2_sparse, prop_h2_sentinel, n_other_sparse) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)
  
  # Sentinel SNP effect
  sentinel_index <- sample(1:n_features, 1)
  beta[sentinel_index] <- rnorm(1, 0, sqrt(h2_sparse * prop_h2_sentinel))
  
  # Remaining heritability for other sparse effects
  h2_sentinel <- h2_sparse * prop_h2_sentinel
  h2_other_sparse <- h2_sparse - h2_sentinel
  
  # Additional sparse SNPs
  if (n_other_sparse > 0) {
    other_sparse_indices <- sample(setdiff(1:n_features, sentinel_index), n_other_sparse)
    beta[other_sparse_indices] <- rnorm(n_other_sparse, 0, sqrt(h2_other_sparse / n_other_sparse))
  } else {
    other_sparse_indices <- integer(0)
  }
  
  # Ensure the sentinel has the largest absolute effect
  if (length(other_sparse_indices) > 0) {
    max_other <- max(abs(beta[other_sparse_indices]))
    if (abs(beta[sentinel_index]) <= max_other) {
      beta[sentinel_index] <- sign(beta[sentinel_index]) * (max_other + 0.01)
    }
  }
  
  # Scale sparse effects to exactly match h2_sparse
  sparse_indices <- c(sentinel_index, other_sparse_indices)
  sparse_effects <- as.vector(G[, sparse_indices] %*% beta[sparse_indices])
  scaling_factor <- sqrt(h2_sparse / var(sparse_effects))
  beta[sparse_indices] <- beta[sparse_indices] * scaling_factor
  
  # Re-check sentinel dominance after scaling
  if (length(other_sparse_indices) > 0) {
    max_other <- max(abs(beta[other_sparse_indices]))
    if (abs(beta[sentinel_index]) <= max_other) {
      beta[sentinel_index] <- sign(beta[sentinel_index]) * (max_other + 0.01)
    }
  }
  
  return(list(beta = beta,
              sentinel_index = sentinel_index,
              other_sparse_indices = other_sparse_indices))
}

###############################################################################
#' Simulate Oligogenic Effects for eQTL Data
#'
#' This function simulates oligogenic effects using a two-component mixture
#' model for a specified number of SNPs not included in the sparse component.
#'
#' @param G A standardized genotype matrix.
#' @param h2_oligogenic Heritability allocated to oligogenic effects.
#' @param n_oligogenic Number of oligogenic SNPs to simulate.
#' @param mixture_props A vector of mixture proportions (must sum to 1).
#' @param mixture_sds A vector of standard deviations for the mixture components.
#' @param non_sparse_indices SNP indices not used in the sparse component.
#' @return A list containing:
#'   \item{beta}{A vector of effect sizes for oligogenic effects (zeros elsewhere).}
#'   \item{oligogenic_indices}{Indices of the oligogenic SNPs.}
#'   \item{mixture_assignments}{A vector (indexed by SNP) of mixture component assignments.}
#' @export
simulate_oligogenic_effects <- function(G, h2_oligogenic, n_oligogenic, mixture_props, mixture_sds, non_sparse_indices) {
  if (abs(sum(mixture_props) - 1) > 1e-6) {
    stop("mixture_props must sum to 1.")
  }
  
  n_features <- ncol(G)
  beta <- rep(0, n_features)
  
  # Ensure we do not select more SNPs than available
  n_available <- length(non_sparse_indices)
  n_oligogenic <- min(n_oligogenic, n_available)
  oligogenic_indices <- sample(non_sparse_indices, n_oligogenic, replace = FALSE)
  
  # Assign mixture components to each oligogenic SNP
  mixture_assignments <- sample(1:length(mixture_props), n_oligogenic, replace = TRUE, prob = mixture_props)
  beta[oligogenic_indices] <- rnorm(n_oligogenic, 0, mixture_sds[mixture_assignments])
  
  # Scale oligogenic effects to match h2_oligogenic
  oligogenic_effects <- as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])
  scaling_factor <- sqrt(h2_oligogenic / var(oligogenic_effects))
  beta[oligogenic_indices] <- beta[oligogenic_indices] * scaling_factor
  
  mixture_assignments_full <- rep(NA, n_features)
  mixture_assignments_full[oligogenic_indices] <- mixture_assignments
  
  return(list(beta = beta,
              oligogenic_indices = oligogenic_indices,
              mixture_assignments = mixture_assignments_full))
}

###############################################################################
#' Simulate Infinitesimal Effects for eQTL Data
#'
#' This function simulates the infinitesimal (polygenic) background effects for
#' a set of SNPs.
#'
#' @param G A standardized genotype matrix.
#' @param h2_infinitesimal Heritability allocated to infinitesimal effects.
#' @param infinitesimal_indices Indices of SNPs to receive infinitesimal effects.
#' @return A vector of effect sizes (with zeros outside infinitesimal_indices).
#' @export
simulate_infinitesimal_effects <- function(G, h2_infinitesimal, infinitesimal_indices) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)
  n_inf <- length(infinitesimal_indices)
  
  if (n_inf > 0) {
    beta[infinitesimal_indices] <- rnorm(n_inf, 0, sqrt(h2_infinitesimal / n_inf))
  }
  return(beta)
}

###############################################################################
#' Generate eQTL Data with Multiple Genetic Architecture Components
#'
#' This function generates simulated gene expression data based on a
#' partitioned genetic architecture model. The total heritability is split into:
#'   - A sparse component (including a dominant sentinel SNP and additional sparse SNPs),
#'   - An oligogenic component (modeled via a two-component mixture),
#'   - An infinitesimal polygenic background.
#'
#' @param G Genotype matrix.
#' @param h2_total Total heritability.
#' @param prop_h2_sparse Proportion of h2_total explained by sparse effects.
#' @param prop_h2_oligogenic Proportion of h2_total explained by oligogenic effects.
#' @param prop_h2_infinitesimal Proportion of h2_total explained by infinitesimal effects.
#' @param prop_h2_sentinel Proportion of h2_sparse explained by the sentinel SNP.
#' @param n_oligogenic Number of oligogenic SNPs to simulate.
#' @param mixture_props Mixture proportions for oligogenic effects (must sum to 1).
#' @param mixture_sds Standard deviations for the mixture components.
#' @param n_other_sparse Number of additional sparse SNPs (besides the sentinel).
#' @param standardize Logical; if TRUE, the genotype matrix will be standardized.
#' @param seed Optional seed for reproducibility.
#' @return A list containing the standardized genotype matrix, simulated phenotype,
#'   combined beta values, indices for each effect component, and realized heritability estimates.
#' @export
generate_eqtl_data <- function(G,
                               h2_total = 0.3,
                               prop_h2_sparse = 0.65,
                               prop_h2_oligogenic = 0.20,
                               prop_h2_infinitesimal = 0.15,
                               prop_h2_sentinel = 0.7,
                               n_oligogenic = 20,
                               mixture_props = c(0.75, 0.25),
                               mixture_sds = c(0.0025, 0.005),
                               n_other_sparse = 2,
                               standardize = TRUE,
                               seed = NULL) {
  # Minimal input validation for proportions
  if (abs(prop_h2_sparse + prop_h2_oligogenic + prop_h2_infinitesimal - 1) > 1e-6) {
    stop("The sum of prop_h2_sparse, prop_h2_oligogenic, and prop_h2_infinitesimal must equal 1.")
  }
  if (abs(sum(mixture_props) - 1) > 1e-6) {
    stop("mixture_props must sum to 1.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Standardize genotype matrix if requested
  if (standardize) {
    G <- scale(G)
  }
  
  n_samples <- nrow(G)
  n_features <- ncol(G)
  
  # Allocate heritability components
  h2_sparse <- h2_total * prop_h2_sparse
  h2_oligogenic <- h2_total * prop_h2_oligogenic
  h2_infinitesimal <- h2_total * prop_h2_infinitesimal
  
  # 1. Sparse Effects (sentinel + additional sparse SNPs)
  sparse_res <- simulate_sparse_effects(G, h2_sparse, prop_h2_sentinel, n_other_sparse)
  beta_sparse <- sparse_res$beta
  sentinel_index <- sparse_res$sentinel_index
  other_sparse_indices <- sparse_res$other_sparse_indices
  sparse_indices <- c(sentinel_index, other_sparse_indices)
  
  # 2. Oligogenic Effects
  non_sparse_indices <- setdiff(1:n_features, sparse_indices)
  oligo_res <- simulate_oligogenic_effects(G, h2_oligogenic, n_oligogenic, mixture_props, mixture_sds, non_sparse_indices)
  beta_oligo <- oligo_res$beta
  oligogenic_indices <- oligo_res$oligogenic_indices
  
  # 3. Infinitesimal Effects (remaining SNPs)
  infinitesimal_indices <- setdiff(non_sparse_indices, oligogenic_indices)
  beta_inf <- simulate_infinitesimal_effects(G, h2_infinitesimal, infinitesimal_indices)
  
  # Combine all effect components
  beta <- beta_sparse + beta_oligo + beta_inf
  
  # Generate latent genetic component and phenotype
  g <- as.vector(G %*% beta)
  var_g <- var(g)
  var_epsilon <- var_g * (1 - h2_total) / h2_total
  epsilon <- rnorm(n_samples, 0, sqrt(var_epsilon))
  y <- g + epsilon
  
  # Calculate realized heritability components
  var_y <- var(y)
  h2_sentinel_actual <- var(G[, sentinel_index] * beta[sentinel_index]) / var_y
  h2_sparse_actual <- var(as.vector(G[, sparse_indices] %*% beta[sparse_indices])) / var_y
  h2_oligogenic_actual <- var(as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])) / var_y
  h2_infinitesimal_actual <- var(as.vector(G[, infinitesimal_indices] %*% beta_inf[infinitesimal_indices])) / var_y
  h2_total_actual <- var(as.vector(G %*% beta)) / var_y
  
  return(list(
    G = G,                # Standardized genotype matrix
    y = y,                # Simulated phenotype (on its generated scale)
    beta = beta,
    h2_total = h2_total_actual,
    h2_sparse = h2_sparse_actual,
    h2_sentinel = h2_sentinel_actual,
    h2_oligogenic = h2_oligogenic_actual,
    h2_infinitesimal = h2_infinitesimal_actual,
    sentinel_index = sentinel_index,
    other_sparse_indices = other_sparse_indices,
    oligogenic_indices = oligogenic_indices,
    infinitesimal_indices = infinitesimal_indices,
    residual_variance = var_epsilon
  ))
}
