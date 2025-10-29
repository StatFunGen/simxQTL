# Check if a matrix is positive definite
#' @keywords internal
is.positive.definite <- function(LD) {
    ev <- eigen(LD, symmetric = TRUE)$values
    return(all(ev > 0))
}

# Get positive definite LD matrix
#' @keywords internal
get_ld_pd <- function(LD, lambda = 1e-3){
    
    if (!isSymmetric(LD))
        stop("Provided LD matrix should be symmetric!")
    
    while (!is.positive.definite(LD)) {
        LD <- LD + lambda * diag(nrow(LD))
        lambda <- lambda * 10
    }
    return(LD)
}

# Calculate pseudo inverse of LD matrix
#' @importFrom Matrix chol2inv
#' @keywords internal
get_ld_inverse <- function(LD, method = "svd", 
                           tol.exact = sqrt(.Machine$double.eps),
                           tol.pct = 1e-3){
    
    if (!isSymmetric(LD))
        stop("Provided LD matrix should be symmetric!")

    if (method == "svd") {
        LD_inv <- get_ld_inv_svd(LD, tol = tol.exact)
    } else if (method == "eigen") {
        LD_inv <- get_ld_inv_eigen(LD, tol = tol.pct)
    } else if (method == "fastchol") {
        LD_inv <- chol2inv(LD)
    } else {
        stop("Invalid method. Choose 'svd', 'eigen', or 'fastchol'.")
    }

    return(LD_inv)

}

# Calculate LD inverse using eigendecomposition
#' @keywords internal
get_ld_inv_eigen <- function(LD, tol = 1e-3){
    
    eigen_LD <- eigen(LD)
    L <- which(cumsum(eigen_LD$values) / sum(eigen_LD$values) > 1-tol)[1]
    LD_inv <- eigen_LD$vectors[,1:L] %*% 
        diag(1/eigen_LD$values[1:L]) %*% 
        t(eigen_LD$vectors[,1:L])
    
    return(LD_inv)
}

# Calculate LD inverse using SVD
#' @keywords internal
get_ld_inv_svd <- function(LD, tol = sqrt(.Machine$double.eps)){
    
    svd_LD <- svd(LD)
    Positive <- svd_LD$d > max(tol * svd_LD$d[1L], 0)
    if (all(Positive))
        LD_inv <- svd_LD$v %*% (1/svd_LD$d * t(svd_LD$u))
    else
        LD_inv <- svd_LD$v[,Positive,drop = FALSE] %*%
        ((1/svd_LD$d[Positive]) * t(svd_LD$u[,Positive,drop = FALSE]))
    
    return(LD_inv)
}
