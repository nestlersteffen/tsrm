
#---- this function generates a permutation matrix 

make_permutationmatrix <- function( n_units=NULL, block_size=NULL, nv=NULL ) 
{
	#- overall dimensionality:
	n_total <- n_units*block_size*nv
	if (nv == 1) return( diag( 1, n_total ) )
	#- make vector with entries for permutation matrix:
	perm <- integer(n_total)
	for (u in 1:n_units) {
		for (v in 1:nv) {
			for (b in 1:block_size) {
	        	from <- (u-1)*nv*block_size + (v-1)*block_size + b  # "Einheit-zuerst"
	        	to   <- (v-1)*n_units*block_size + (u-1)*block_size + b  # "Variable-zuerst"
	        	perm[to] <- from
	      	}
	    }
	}
	P <- matrix(0, n_total, n_total)
	for (j in 1:n_total) P[j, perm[j]] <- 1
	return( P )
}
