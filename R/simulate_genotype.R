#' Simulate Independent Genotypes
#'
#' @param n Sample size.
#' @param p Number of SNPs.
#' @param min_maf Minimum minor allele frequency.
#' @param max_maf Maximum minor allele frequency.
#' @param scale Logical, whether to scale the data.
#' @return A matrix of genotypes.
#' @examples
#' sim_geno_indep(n = 100, p = 10, min_maf = 0.01, max_maf = 0.4, scale = TRUE)
#' @importFrom purrr map2_dfc
#' @export
sim_geno_indep <- function(n, p, min_maf = 0.01, max_maf = 0.4, scale = FALSE) {
  if (missing(n)) stop("Please provide the sample size")
  if (is.null(p)) stop("Please provide the number of independent SNPs!")

  mafs <- runif(p, min_maf, max_maf)
  G <- map2_dfc(.x = mafs, .y = rep(n, p), ~ rbinom(.y, size = 2, prob = .x))
  colnames(G) <- paste0("variant", 1:p)

  if (scale) {
    G <- scale(G)
  }

  return(G)
}

#' Simulate Genotypes Based on LD Matrix
#'
#' @param n Sample size.
#' @param LD LD matrix.
#' @param min_maf Minimum minor allele frequency.
#' @param max_maf Maximum minor allele frequency.
#' @param lambda Regularization parameter.
#' @param is.discrete Logical, whether to generate discrete (TRUE) or continuous (FALSE) genotypes.
#' @param scale Logical, whether to scale the data.
#' @param tol Tolerance for checking diagonal elements of LD matrix.
#' @return A matrix of genotypes.
#' @examples
#' sim_geno_LD(n = 100, LD = matrix(runif(100), 10, 10), min_maf = 0.01, max_maf = 0.4, lambda = 1e-3, is.discrete = FALSE, scale = TRUE)
#' @export
sim_geno_LD <- function(n, LD, min_maf = 0.01, max_maf = 0.4, lambda = 1e-3, is.discrete = FALSE, scale = FALSE, tol = sqrt(.Machine$double.eps)) {
  if (missing(n)) stop("Please provide the sample size")
  if (is.null(LD)) stop("Please provide LD matrix!")

  if (!isSymmetric(LD)) stop("LD matrix should be symmetric!")
  if (any(diag(LD) < 1 - tol | diag(LD) > 1 + tol)) stop("Diagonal elements of LD matrix should be 1!")

  LD <- as.matrix(LD)
  p <- nrow(LD)
  maf <- runif(p, min = min_maf, max = max_maf)

  G <- safe_rmvnorm(n, p, LD, maf, lambda)
  if (is.discrete) {
      G <- binomialize_genotype(G, p, maf)
  }

  if (scale) {
    G <- scale(G)
  }

  colnames(G) <- paste0("variant", 1:p)
  return(G)
}

#' Simulate Genotypes Based on Real Data
#'
#' @param n Sample size.
#' @param file_path Path to the genotype file.
#' @param min_maf Minimum minor allele frequency.
#' @param scale Logical, whether to scale the data.
#' @return A matrix of genotypes.
#' @examples
#' sim_geno_real(n = 100, file_path = "path/to/real/file", min_maf = 0.01, scale = TRUE)
#' @export
sim_geno_real <- function(n, file_path, min_maf = 0.01, scale = FALSE) {
  if (missing(n)) stop("Please provide the sample size")
  if (is.null(file_path)) stop("Please provide the path to genotype data!")
  if (!grepl("\\.bed$", file_path)) stop("Please provide plink bfiles!")

  G <- process_ukb(file_path, n, min_maf)

  if (scale) {
    G <- scale(G)
  }

  colnames(G) <- paste0("variant", 1:ncol(G))
  return(G)
}

# Safe multivariate normal random generation
#' @importFrom mvtnorm rmvnorm
#' @keywords internal
safe_rmvnorm <- function(n, p, LD, maf, lambda = 1e-3) {
    
    var_g <- diag(sqrt(2*maf*(1-maf)))
    Sigma <- var_g %*% LD %*% var_g
    G <- try(rmvnorm(n, mean = rep(0, p), sigma = Sigma), silent = TRUE)
    if(any(class(G) == "try-error")) {
        Sigma <- get_ld_pd(Sigma, lambda = lambda)
        G <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    }
    return(G)
    
}

# Convert continuous genotypes to discrete (binomial)
#' @keywords internal
binomialize_genotype <- function(G, p, maf) {
    
    sapply(1:p, function(i) {
        tmp.maf <- maf[i]
        g <- G[,i]
        tmp <- sort(g, decreasing = TRUE)
        p2 <- ceiling(p * tmp.maf^2)
        p1 <- ceiling(p - p * (1 - tmp.maf)^2)
        2 * (g >= tmp[p2]) + 1 * (g >= tmp[p1 + p2]) * (g < tmp[p2])
    })
    
}

# Process UK Biobank genotype file
#' @importFrom BEDMatrix BEDMatrix
#' @importFrom data.table fread
#' @keywords internal
process_ukb <- function(file.path, n, min_maf = 0.01) {
    
    G <- BEDMatrix(file.path)
    G <- data.matrix(G[1:n,])
    G <- apply(G, 2, function(g) {
        tmp <- which(is.na(g))
        g[tmp] <- mean(g, na.rm = TRUE)
        g
    })
    maf <- colMeans(G, na.rm = TRUE) / 2
    snpname <- strsplit(file.path, split = "\\.bed")[[1]][1]
    snp.names <- unlist(fread(paste0(snpname, ".bim"))[,2])
    colnames(G) <- snp.names
    poss <- which(maf <= min_maf)
    if (length(poss) != 0) {
        G <- G[, -poss, drop = FALSE]
    }
    return(G)
    
}


# Generate SNP with specified LD
#' @keywords internal
generate_snp_LD <- function(g, ld, maf) {
    
    # - change genotype to haplotype
    haplo_transform <- function(i.g) {
        switch(as.character(i.g),
               "0" = c(0, 0),
               "2" = c(1, 1),
               {
                   g1 <- as.integer(runif(1) > 0.5)
                   c(g1, 1 - g1)
               })
    }
    g.hap <- t(sapply(g, haplo_transform))
    
    # - generate haplotype based on maf and ld
    n <- length(g)
    g2 <- matrix(rbinom(2 * n, 1, maf), nrow = n)
    
    update_g2 <- function(g.hap.col, g2.col) {
        pA <- mean(g.hap.col)
        pAB <- pA * maf + sqrt(ld^2 * pA * (1 - pA) * maf * (1 - maf))
        pBA <- pAB / pA
        pbA <- (pA * maf - pAB) / pA
        pBa <- (maf - pAB) / (1 - pA)
        pba <- (1 - pA - maf + pAB) / (1 - pA)
        
        sample(c(0, 1), n, replace = TRUE, prob = ifelse(g.hap.col == 1, c(pbA, pBA), c(pba, pBa)))
    }
    
    g2[, 1] <- update_g2(g.hap[, 1], g2[, 1])
    g2[, 2] <- update_g2(g.hap[, 2], g2[, 2])
    
    rowSums(g2)
}
