###############################################################################
#' Simulate Sparse Effects for eQTL Data
#'
#' This function simulates the sparse component of the genetic effects using
#' target-based scaling. Effects are drawn from a normal distribution and then
#' scaled to match the target heritability.
#'
#' @param G A standardized genotype matrix.
#' @param h2_sparse Heritability allocated to sparse effects.
#' @param n_sparse Number of sparse SNPs.
#' @param effect_sd Standard deviation for drawing sparse effects (default 0.5).
#' @return A list containing:
#'   \item{beta}{A vector of effect sizes with nonzero entries for the sparse SNPs.}
#'   \item{sparse_indices}{Indices of the sparse SNPs.}
#' @keywords internal
simulate_sparse_effects <- function(G, h2_sparse, n_sparse, effect_sd = 0.5) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)

  # Select sparse SNPs
  if (n_sparse > 0) {
    sparse_indices <- sample(1:n_features, n_sparse)

    # Draw effects from normal distribution with specified SD
    beta[sparse_indices] <- rnorm(n_sparse, mean = 0, sd = effect_sd)
  } else {
    sparse_indices <- integer(0)
  }

  # Scale sparse effects to exactly match h2_sparse
  if (length(sparse_indices) > 0) {
    sparse_effects <- as.vector(G[, sparse_indices, drop = FALSE] %*% beta[sparse_indices])
    current_var <- var(sparse_effects)

    if (current_var > 0) {
      scaling_factor <- sqrt(h2_sparse / current_var)
      beta[sparse_indices] <- beta[sparse_indices] * scaling_factor
    }
  }

  return(list(beta = beta,
              sparse_indices = sparse_indices))
}

###############################################################################
#' Simulate Oligogenic Effects for eQTL Data
#'
#' This function simulates oligogenic effects using target-based scaling with
#' a two-component mixture model. Effects are drawn from a mixture distribution
#' and then scaled to match the target heritability.
#'
#' @param G A standardized genotype matrix.
#' @param h2_oligogenic Heritability allocated to oligogenic effects.
#' @param n_oligogenic Number of oligogenic SNPs to simulate.
#' @param mixture_props A vector of mixture proportions (must sum to 1).
#' @param non_sparse_indices SNP indices not used in the sparse component.
#' @param effect_sds Standard deviations for the two mixture components (default c(0.05, 0.15)).
#' @return A list containing:
#'   \item{beta}{A vector of effect sizes for oligogenic effects (zeros elsewhere).}
#'   \item{oligogenic_indices}{Indices of the oligogenic SNPs.}
#'   \item{mixture_assignments}{A vector (indexed by SNP) of mixture component assignments.}
#' @keywords internal
simulate_oligogenic_effects <- function(G, h2_oligogenic, n_oligogenic, mixture_props,
                                        non_sparse_indices,
                                        effect_sds = c(0.05, 0.15)) {
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

  # Draw effects from mixture components
  beta[oligogenic_indices] <- rnorm(n_oligogenic, mean = 0, sd = effect_sds[mixture_assignments])

  # Scale oligogenic effects to exactly match h2_oligogenic
  oligogenic_effects <- as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])
  current_var <- var(oligogenic_effects)

  if (current_var > 0) {
    scaling_factor <- sqrt(h2_oligogenic / current_var)
    beta[oligogenic_indices] <- beta[oligogenic_indices] * scaling_factor
  }

  mixture_assignments_full <- rep(NA, n_features)
  mixture_assignments_full[oligogenic_indices] <- mixture_assignments

  return(list(beta = beta,
              oligogenic_indices = oligogenic_indices,
              mixture_assignments = mixture_assignments_full))
}

###############################################################################
#' Simulate Infinitesimal Effects for eQTL Data
#'
#' This function simulates the infinitesimal (polygenic) background effects using
#' target-based scaling. Effects are drawn from a normal distribution and scaled
#' to match the target heritability. The small effect_sd ensures infinitesimal
#' effects remain tiny.
#'
#' @param G A standardized genotype matrix.
#' @param h2_infinitesimal Heritability allocated to infinitesimal effects.
#' @param infinitesimal_indices Indices of SNPs to receive infinitesimal effects.
#' @param effect_sd Standard deviation for drawing infinitesimal effects (default 0.01).
#' @return A vector of effect sizes (with zeros outside infinitesimal_indices).
#' @keywords internal
simulate_infinitesimal_effects <- function(G, h2_infinitesimal, infinitesimal_indices,
                                           effect_sd = 0.01) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)
  n_inf <- length(infinitesimal_indices)

  if (n_inf > 0) {
    # Draw from normal distribution with small SD
    beta[infinitesimal_indices] <- rnorm(n_inf, mean = 0, sd = effect_sd)

    # Scale to exactly match h2_infinitesimal
    inf_effects <- as.vector(G[, infinitesimal_indices] %*% beta[infinitesimal_indices])
    current_var <- var(inf_effects)

    if (current_var > 0) {
      scaling_factor <- sqrt(h2_infinitesimal / current_var)
      beta[infinitesimal_indices] <- beta[infinitesimal_indices] * scaling_factor
    }
  }
  return(beta)
}

###############################################################################
#' Identify Causal SNPs Based on Statistical Power
#'
#' This function identifies SNPs that have sufficient statistical power to be
#' detected as causal. Power is calculated based on the non-centrality parameter
#' of a chi-square test, considering sample size, effect size, and SNP variance.
#'
#' @param G Genotype matrix (samples × SNPs).
#' @param beta Vector of true effect sizes for each SNP.
#' @param residual_variance Residual (error) variance.
#' @param power Minimum power threshold for a SNP to be considered causal (default 0.80).
#' @return A vector of indices corresponding to SNPs with power >= threshold.
#' @details The function calculates power using a non-central chi-square distribution
#'   with non-centrality parameter NCP = n * beta^2 * var(X) / residual_variance, where the
#'   significance threshold is Bonferroni-corrected (alpha = 0.05 / p).
#' @keywords internal
is_causal_power <- function(G, beta, residual_variance, power = 0.80) {
  n <- nrow(G)
  p <- ncol(G)
  alpha <- 0.05 / p

  var_snp <- apply(G, 2, var)
  ncp <- n * (beta^2) * var_snp / residual_variance

  crit_val <- qchisq(alpha, df = 1, lower.tail = FALSE)
  power_per_snp <- pchisq(crit_val, df = 1, ncp = ncp, lower.tail = FALSE)

  causal_idx <- which(power_per_snp >= power)
  return(causal_idx)
}

###############################################################################
#' Generate eQTL Data with Multiple Genetic Architecture Components
#'
#' This function generates simulated gene expression data with a partitioned
#' genetic architecture that enforces strict effect size hierarchies:
#' |sparse| > |oligogenic| >> |infinitesimal|
#'
#' @param G Genotype matrix.
#' @param h2g Total SNP heritability (proportion of variance explained by genotyped SNPs).
#' @param prop_h2_sparse Proportion of h2g explained by sparse effects.
#' @param prop_h2_oligogenic Proportion of h2g explained by oligogenic effects.
#' @param prop_h2_infinitesimal Proportion of h2g explained by infinitesimal effects.
#' @param n_sparse Number of sparse SNPs.
#' @param n_oligogenic Number of oligogenic SNPs to simulate.
#' @param mixture_props Mixture proportions for oligogenic effects (must sum to 1). Default c(0.75, 0.25) means 75% smaller effects, 25% larger effects.
#' @param sparse_sd Standard deviation for drawing sparse effects (default 0.5).
#' @param oligo_sds Standard deviations for oligogenic mixture components (default c(0.05, 0.15)).
#' @param inf_sd Standard deviation for drawing infinitesimal effects (default 0.01).
#' @param standardize Logical; if TRUE, the genotype matrix will be standardized.
#' @param independent Logical; if TRUE, ensures all sparse and oligogenic SNPs have |r| < 0.15 with each other (default FALSE).
#' @param verbose Logical; if TRUE, prints progress messages including LD constraint attempts (default FALSE).
#' @param seed Optional seed for reproducibility.
#' @return A list containing the standardized genotype matrix, simulated phenotype,
#'   combined beta values, indices for each effect component, realized heritability estimates,
#'   effect size ranges, hierarchy validation results, and causal indices.
#' @export
generate_eqtl_data <- function(G,
                               h2g = 0.30,
                               prop_h2_sparse = 0.50,
                               prop_h2_oligogenic = 0.15,
                               prop_h2_infinitesimal = 0.35,
                               n_sparse = 2,
                               n_oligogenic = 10,
                               mixture_props = c(0.75, 0.25),
                               sparse_sd = 0.5,
                               oligo_sds = c(0.05, 0.15),
                               inf_sd = 0.01,
                               standardize = TRUE,
                               independent = TRUE,
                               verbose = FALSE,
                               seed = NULL) {
  # Input validation
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
  h2_sparse <- h2g * prop_h2_sparse
  h2_oligogenic <- h2g * prop_h2_oligogenic
  h2_infinitesimal <- h2g * prop_h2_infinitesimal

  # Setup for LD constraint checking
  max_attempts <- if (independent) 200 else 1
  attempt <- 0
  ld_satisfied <- FALSE

  # Main data generation loop
  if (verbose && independent) {
    cat("Attempting to satisfy LD constraints:\n")
    cat("  - Sparse SNPs: |r| < 0.15\n")
    cat("  - Oligogenic SNPs: |r| < 0.15\n")
    cat("  - Sparse-Oligogenic cross: |r| < 0.15\n")
    cat("  (Safeguard: will accept results after 200 attempts)\n")
  }

  while (!ld_satisfied && attempt < max_attempts) {
    attempt <- attempt + 1

    if (verbose && independent && attempt > 1) {
      cat("  Attempt", attempt, "of", max_attempts, "...\n")
    }

    # Set seed for this attempt (for reproducibility)
    if (!is.null(seed) && attempt > 1) {
      set.seed(seed * 1000 + attempt)
    }

    # 1. Sparse Effects using target-based scaling
    sparse_res <- simulate_sparse_effects(G, h2_sparse, n_sparse, sparse_sd)
    beta_sparse <- sparse_res$beta
    sparse_indices <- sparse_res$sparse_indices

    # 2. Oligogenic Effects using target-based scaling
    non_sparse_indices <- setdiff(1:n_features, sparse_indices)
    oligo_res <- simulate_oligogenic_effects(G, h2_oligogenic, n_oligogenic,
                                             mixture_props, non_sparse_indices,
                                             oligo_sds)
    beta_oligo <- oligo_res$beta
    oligogenic_indices <- oligo_res$oligogenic_indices

    # 3. Infinitesimal Effects using target-based scaling
    infinitesimal_indices <- setdiff(non_sparse_indices, oligogenic_indices)
    beta_inf <- simulate_infinitesimal_effects(G, h2_infinitesimal, infinitesimal_indices,
                                               inf_sd)

    # Combine all effect components
    beta <- beta_sparse + beta_oligo + beta_inf

    # Generate latent genetic component and phenotype
    g <- as.vector(G %*% beta)
    var_g <- var(g)
    var_epsilon <- var_g * (1 - h2g) / h2g
    epsilon <- rnorm(n_samples, 0, sqrt(var_epsilon))
    y <- g + epsilon

    # Check LD constraint if independent = TRUE (applied to sparse and oligogenic SNPs)
    if (independent) {
      ld_satisfied <- TRUE  # Assume satisfied unless violations found

      # Check sparse SNPs: |r| < 0.15
      if (length(sparse_indices) > 1) {
        sparse_cor <- cor(G[, sparse_indices, drop = FALSE])
        high_ld_sparse <- which(abs(sparse_cor) >= 0.15 & upper.tri(sparse_cor, diag = FALSE), arr.ind = TRUE)

        if (nrow(high_ld_sparse) > 0) {
          ld_satisfied <- FALSE
        }
      }

      # Check oligogenic SNPs: |r| < 0.15 (changed from 0.30)
      if (ld_satisfied && length(oligogenic_indices) > 1) {
        oligo_cor <- cor(G[, oligogenic_indices, drop = FALSE])
        high_ld_oligo <- which(abs(oligo_cor) >= 0.15 & upper.tri(oligo_cor, diag = FALSE), arr.ind = TRUE)

        if (nrow(high_ld_oligo) > 0) {
          ld_satisfied <- FALSE
        }
      }

      # Check between sparse and oligogenic: |r| < 0.15
      if (ld_satisfied && length(sparse_indices) > 0 && length(oligogenic_indices) > 0) {
        cross_cor <- cor(G[, sparse_indices, drop = FALSE], G[, oligogenic_indices, drop = FALSE])
        high_ld_cross <- which(abs(cross_cor) >= 0.15, arr.ind = TRUE)

        if (nrow(high_ld_cross) > 0) {
          ld_satisfied <- FALSE
        }
      }
    } else {
      # If not independent, accept
      ld_satisfied <- TRUE
    }
  }

  # Report LD constraint results
  if (independent && !ld_satisfied) {
    msg <- paste0("Failed to satisfy LD constraints after ",
                  max_attempts, " attempts. Returning last generated data with LD violations.")
    warning(msg)
    if (verbose) {
      cat("\n", msg, "\n")
      # Report actual LD values
      if (length(sparse_indices) > 1) {
        sparse_cor <- cor(G[, sparse_indices, drop = FALSE])
        max_sparse_cor <- max(abs(sparse_cor[upper.tri(sparse_cor)]))
        cat("  Final max |r| among sparse SNPs:", round(max_sparse_cor, 4),
            ifelse(max_sparse_cor < 0.15, "✓", "✗ (constraint: < 0.15)"), "\n")
      }
      if (length(oligogenic_indices) > 1) {
        oligo_cor <- cor(G[, oligogenic_indices, drop = FALSE])
        max_oligo_cor <- max(abs(oligo_cor[upper.tri(oligo_cor)]))
        cat("  Final max |r| among oligogenic SNPs:", round(max_oligo_cor, 4),
            ifelse(max_oligo_cor < 0.15, "✓", "✗ (constraint: < 0.15)"), "\n")
      }
      if (length(sparse_indices) > 0 && length(oligogenic_indices) > 0) {
        cross_cor <- cor(G[, sparse_indices, drop = FALSE], G[, oligogenic_indices, drop = FALSE])
        max_cross_cor <- max(abs(cross_cor))
        cat("  Final max |r| between sparse and oligogenic:", round(max_cross_cor, 4),
            ifelse(max_cross_cor < 0.15, "✓", "✗ (constraint: < 0.15)"), "\n")
      }
    }
  } else if (verbose && independent) {
    if (attempt == 1) {
      cat("✓ LD constraints satisfied on first attempt!\n")
    } else {
      cat("✓ LD constraints satisfied after", attempt, "attempts.\n")
    }
  }

  # Calculate causal indices using power-based identification
  causal_indices <- is_causal_power(G, beta, var_epsilon, power = 0.80)

  # Calculate empirical heritability for each component (final data)
  var_y <- var(y)
  h2_sparse_actual <- var(as.vector(G[, sparse_indices] %*% beta[sparse_indices])) / var_y
  h2_oligogenic_actual <- var(as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])) / var_y
  h2_infinitesimal_actual <- var(as.vector(G[, infinitesimal_indices] %*% beta_inf[infinitesimal_indices])) / var_y
  h2g_actual <- var(as.vector(G %*% beta)) / var_y

  return(list(
    G = G,
    y.ori = y,
    y = scale(y),
    beta = beta,
    h2g = h2g_actual,
    h2_sparse = h2_sparse_actual,
    h2_oligogenic = h2_oligogenic_actual,
    h2_infinitesimal = h2_infinitesimal_actual,
    sparse_indices = sparse_indices,
    oligogenic_indices = oligogenic_indices,
    infinitesimal_indices = infinitesimal_indices,
    residual_variance = var_epsilon,
    causal = causal_indices
  ))
}
