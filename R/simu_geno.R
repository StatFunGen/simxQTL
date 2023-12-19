
#' Simulate Genotypes
#'
#' Sample genotypes from a multivariate normal approximation.
#'
# method = "LD": generate genotype matrix based on LD matrix
# method = "UKB": generate genotype matrix based on sim1000G package
# method = "indep": generate independent genotype matrix

# --- potential issue: 0,1,2 genotype can not keep LD.

sim_geno <- function(n, file.path = NULL, p = NULL, LD = NULL, 
                     method = "LD", type = "continous",
                     min_maf = 0.01, max_maf = 0.4,
                     tol = 1e-9, lambda = 1e-3) {
    
    if (missing(n)) stop("Please provide the sample size")
    
    switch(method,
           indep = {
               if (is.null(p)) 
                   stop("Please provide the number of independent SNPs!")
               G <- matrix(rbinom(n * p, 2, runif(n * p, min_maf, max_maf)), ncol = p)
               colnames(G) <- paste0("variant", 1:p)
           },
           LD = {
               if (is.null(file.path) && is.null(LD)) 
                   stop("Please provide LD matrix or path of LD structure!")
               
               if (!is.null(file.path)) {
                   LD <- readRDS(file.path)
               } else {
                   if (!isSymmetric(LD)) stop("LD matrix should be symmetric!")
                   if (any(diag(LD) < 1 - tol | diag(LD) > 1 + tol))
                       stop("Diagonal elements of LD matrix should be 1!")
               }
               
               LD <- as.matrix(LD)
               if (is.null(colnames(LD))) {
                   snp.names <- paste0("variant", 1:nrow(LD))
               } else {
                   snp.names <- colnames(LD)
               }
               p <- nrow(LD)
               maf <- runif(p, min = min_maf, max = max_maf)
               
               G <- safe_rmvnorm(n, p, LD, maf, lambda)
               if (type == "discrete") {
                   G <- binomialize_genotype(G, p, maf)
                   # ---- if use the haplotype procedure 
                   #     G <- matrix(0, nrow = n, ncol = p)
                   #     G[,1] <- rbinom(n, 2, maf[1])
                       # for (i in 2:p){
                       #     g <- G[,i-1]
                       #     G[,i] <- generate_snp_LD(g, LD[i-1, i], maf[i])
                       # }
               }
               colnames(G) <- snp.names
           },
           UKB = {
               if (is.null(file.path)) stop("Please provide the path of UK Biobank!")
               if (!grepl("\\.bed$", file.path)) stop("Please provide plink bfiles!")
               G <- process_ukb(file.path, n, min_maf)
           },
           stop("Invalid method! Valid methods are 'indep', 'LD', 'UKB'")
    )
    
    return(G)
}

safe_rmvnorm <- function(n, p, LD, maf, lambda = 1e-3) {
    
    var_g <- diag(sqrt(2*maf*(1-maf)))
    Sigma <- var_g %*% LD %*% var_g
    G <- try(mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma), silent = TRUE)
    if(any(class(G) == "try-error")) {
        Sigma <- get_ld_pd(Sigma, lambda = lambda)
        G <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    }
    return(G)
    
}

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

process_ukb <- function(file.path, n, min_maf = 0.01) {
    
    G <- BEDMatrix::BEDMatrix(file.path)
    G <- data.matrix(G[1:n,])
    G <- apply(G, 2, function(g) {
        tmp <- which(is.na(g))
        g[tmp] <- mean(g, na.rm = TRUE)
        g
    })
    maf <- colMeans(G, na.rm = TRUE) / 2
    snpname <- strsplit(file.path, split = "\\.bed")[[1]][1]
    snp.names <- unlist(data.table::fread(paste0(snpname, ".bim"))[,2])
    colnames(G) <- snp.names
    poss <- which(maf <= min_maf)
    if (length(poss) != 0) {
        G <- G[, -poss, drop = FALSE]
    }
    return(G)
    
}


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
