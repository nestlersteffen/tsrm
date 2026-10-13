

bayes_dhalfcauchy <- function( x, s = 2.5, log = FALSE ) {
	# dens <- log(2 * s) - log(pi * {
    #     x * x + s * s
    # })
    dens <- log(2) - log(pi*s) - log(1 + (x/s)^2)
	if ( log == FALSE )
		dens <- exp( dens )
	return( dens )
}

bayes_dhalft <- function (x, scale = 25, nu = 1, log = FALSE) 
{
    x <- as.vector(x)
    scale <- as.vector(scale)
    nu <- as.vector(nu)
    if (any(scale <= 0)) 
        stop("The scale parameter must be positive.")
    NN <- max(length(x), length(scale), length(nu))
    x <- rep(x, len = NN)
    scale <- rep(scale, len = NN)
    nu <- rep(nu, len = NN)
    dens <- log(2) - log(scale) + lgamma((nu + 1)/2) - lgamma(nu/2) - 
        0.5 * log(pi * nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/scale) * 
        (x/scale))
    if (log == FALSE) 
        dens <- exp(dens)
    return(dens)
}

# lkj-density omits normalization constant
bayes_dlkjcorr <- function( R , eta = 1 , log=TRUE ) {
    dens <- det(R)^(eta-1)
    if ( log == TRUE ) dens <- log( dens )
    return( dens )
}

# function to generate a draw from LKJ distribution
# taken from Ben Bolker, K is dimension of matrix
bayes_rlkjcorr <- function ( n , K , eta = 1 ) {

    stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
    stopifnot(eta > 0)
    #if (K == 1) return(matrix(1, 1, 1))
    
    f <- function() {
        alpha <- eta + (K - 2)/2
        r12 <- 2 * rbeta(1, alpha, alpha) - 1
        R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
        R[1,1] <- 1
        R[1,2] <- r12
        R[2,2] <- sqrt(1 - r12^2)
        if(K > 2) for (m in 2:(K - 1)) {
            alpha <- alpha - 0.5
            y <- rbeta(1, m / 2, alpha)
            
            # Draw uniformally on a hypersphere
            z <- rnorm(m, 0, 1)
            z <- z / sqrt(crossprod(z)[1])
            
            R[1:m,m+1] <- sqrt(y) * z
            R[m+1,m+1] <- sqrt(1 - y)
        }
        return(crossprod(R))
    }
    R <- replicate( n , f() )
    if ( dim(R)[3]==1 ) {
        R <- R[,,1]
    } else {
        # need to move 3rd dimension to front, so conforms to array structure that Stan uses
        R <- aperm(R,c(3,1,2))
    }
    return(R)
}