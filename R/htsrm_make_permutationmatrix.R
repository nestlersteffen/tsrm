
#---- this function generates a permutation matrix

htsrm_make_permutationmatrix <- function( n_units=NULL, block_size=NULL, nv=NULL )
{
	#- scalar -> replicate to length nv (backward compatible):
	if ( length( block_size ) == 1L ) {
		block_size <- rep( block_size, nv )
	}
	if ( length( block_size ) != nv ) {
		stop( "block_size must be a scalar or a vector of length nv." )
	}

	#- overall dimensionality:
	S       <- sum( block_size )      # effects per unit, summed over variables
	n_total <- n_units * S
	if ( nv == 1 ) return( diag( 1, n_total ) )

	#- cumulative offsets (replace the constant (v-1)*block_size terms):
	cs_within <- c( 0, cumsum( block_size ) )            # offset WITHIN a unit (unit-major)
	cs_var    <- c( 0, cumsum( n_units * block_size ) )  # offset of the variable block (var-major)

	#- make vector with entries for permutation matrix:
	perm <- perm_inv <- integer( n_total )
	for ( u in 1:n_units ) {
		for ( v in 1:nv ) {
			s_v <- block_size[v]
			if ( s_v == 0 ) next                                  # variable absent on this level
			for ( b in 1:s_v ) {
				from <- (u-1)*S + cs_within[v] + b                # "Einheit-zuerst" (unit-major)
				to   <- cs_var[v] + (u-1)*s_v + b                 # "Variable-zuerst" (var-major)
				perm[to]       <- from
				perm_inv[from] <- to
			}
		}
	}
	P <- matrix( 0, n_total, n_total )
	for ( j in 1:n_total ) P[j, perm[j]] <- 1
	return( list( perm = perm, perm_inv = perm_inv, P = P ) )
}