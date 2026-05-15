

srm_make_bldiag <- function( mat_list = NULL ) {
	
	#- get relevant dimnesions:
	ncolsum <- Reduce("+", lapply( mat_list, function(x) ncol(x) ) )
	nrowsum <- Reduce("+", lapply( mat_list, function(x) nrow(x) ) )

	#- make block matrix:
	block <- matrix(0, nrow=nrowsum, ncol=ncolsum )

	#- fill in the matrices:
	idx_row <- 0
	idx_col <- 0
	
	for ( ii in 1:length( mat_list ) ) {

		#- get matrix in list:
		mat <- mat_list[[ii]]
		nrowmat <- nrow(mat)
		ncolmat <- ncol(mat)

		#- insert
		block[(idx_row+1):(idx_row+nrowmat),(idx_col+1):(idx_col+ncolmat)] <- mat


		#- increase idx-variables:
		idx_row <- idx_row + nrowmat 
		idx_col <- idx_col + ncolmat

	}

	#- return result
	return( block )
}