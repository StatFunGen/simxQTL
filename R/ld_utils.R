

# - Calculate LD matrix
get_ld <- function(X, intercepte = FALSE){
    X = t(X)
    # Center each variable
    if (!intercepte){
        X = X - rowMeans(X)
    }
    # Standardize each variable
    X = X / sqrt(rowSums(X^2))
    # Calculate correlations
    cr = tcrossprod(X)
    return(cr)
}

# - Get positive definite LD matrix
# Function to check if a matrix is positive definite
is.positive.definite <- function(LD) {
    ev <- eigen(LD, symmetric = TRUE)$values
    return(all(ev > 0))
}

get_ld_pd <- function(LD, lambda = 1e-3){
    
    if (!isSymmetric(LD))
        stop("Provided LD matrix should be sysmatric!")
    
    while (!is.positive.definite(LD)) {
        LD <- LD + lambda * diag(nrow(LD))
        lambda <- lambda * 10
    }
    return(LD)
}



# - Calculate pseduo inverse of LD matrix
get_ld_inverse <- function(LD, method = "svd", 
                           tol.exact = sqrt(.Machine$double.eps),
                           tol.pct = 1e-3){
    
    if (!isSymmetric(LD))
        stop("Provided LD matrix should be sysmatric!")
    
    if (method == "svd")
        LD_inv <- get_ld_inv_svd(LD, tol = tol.exact)
    
    if (method == "eigen")
        LD_inv <- get_ld_inv_eigen(LD, tol = tol.pct)
    
    if (method == "fastchol")
        LD_inv <- Matrix::chol2inv(LD)
    
    return(LD_inv)
    
}

get_ld_inv_eigen <- function(LD, tol = 1e-3){
    
    eigen_LD <- eigen(LD)
    L <- which(cumsum(eigen_LD$values) / sum(eigen_LD$values) > 1-tol)[1]
    LD_inv <- eigen_LD$vectors[,1:L] %*% 
        diag(1/eigen_LD$values[1:L]) %*% 
        t(eigen_LD$vectors[,1:L])
    
    return(LD_inv)
}

get_ld_inv_svd <- function(LD, tol = sqrt(.Machine$double.eps)){
    
    svd_LD <- svd(LD)
    Positive <- svd_LD$d > max(tol * svd_LD$d[1L], 0)
    if (all(Positive))
        LD_inv = svd_LD$v %*% (1/svd_LD$d * t(svd_LD$u))
    else 
        LD_inv = svd_LD$v[,Positive,drop = FALSE] %*%
        ((1/svd_LD$d[Positive]) * t(svd_LD$u[,Positive,drop = FALSE]))
    
    return(LD_inv)
}

