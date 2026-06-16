
make_pd <- function(H, tol_factor = 1e-6) {
    eig        <- eigen(H, symmetric = TRUE)
    eig$values <- pmax(eig$values, tol_factor )
    H_pd       <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    return(H_pd)
}