
#---- function to compute the update of BETA in TSRM-Bayes 

# tsrm_bayes_compute_Vepsi <- function( L_T=NULL, Zti=NULL, p=NULL, nt=NULL)
# {
# 	# define matrices
# 	ni <- nrow(Zti)
# 	B  <- aperm(array(Zti, c(ni, p, nt)), c(1, 3, 2))   # ni x nt x p
# 	dim(B) <- c(ni * nt, p)
# 	B <- B %*% L_T                                      # alle Blöcke Zk %*% L auf einmal
# 	dim(B) <- c(ni, nt, p)
# 	B <- aperm(B, c(1, 3, 2))
# 	dim(B) <- c(ni, p * nt)
# 	return( tcrossprod(B) )
# }

tsrm_bayes_compute_Vepsi <- function( L_T=NULL, Zti=NULL, nt=NULL)
{
	# define matrices
	p  <- ncol( L_T )
	ni <- nrow(Zti)
	B  <- matrix(0, ni, nt * p)
	for (k in seq_len(nt)) {
	  idx <- (k - 1) * p + seq_len(p)
	  B[, idx] <- Zti[, idx, drop = FALSE] %*% L_T
	}
	return( tcrossprod(B) )
}

tsrm_bayes_update_beta <- function( data_list=NULL, args_list=NULL, person=NULL, dyads=NULL, 
	SIGMA_T=NULL ) 
{

	#- get matrices and groupinfo:
  	groupinfo <- data_list[["groupinfo"]]
  	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	Zt <- data_list[["Zt"]]
	nt <- data_list[["nt"]]

	#- get no. groups:
	ngroups <- nrow( groupinfo )

	#- compute G_T:
	if ( !args_list$large ) {
		G_T <- L_T <- diag( 1, nt ) %x% SIGMA_T
	} else {
		G_T <- L_T <- t(chol(SIGMA_T))
	}
	
	#- get hyperparms:
	mb 	  <- args_list$mcmc_args$mb
	invMb <- args_list$mcmc_args$invMb

	if ( args_list$use_rcpp ) {

		#- change groupinfo to integer:
		storage.mode( groupinfo ) <- "integer"

		#- call the C++ function
		tmp <- tsrm_bayes_update_beta_export( y=y, X=X, Zp=Zp, Zd=Zd, Zt=Zt, nt=nt, 
			G_T=G_T, L_T=L_T, groupinfo=groupinfo, ngroups=ngroups, person=person, 
			dyads=dyads, mb=mb, invMb=invMb, args_list$large )
		out <- list( BETA = as.matrix( tmp$BETA ), OBJ_invVepsi = tmp$OBJ_invVepsi )

	} else {

		#- initialize:
		sum_mu <- matrix( 0, nrow = ncol( X ), ncol = 1 )
		sum_om <- diag( 0, ncol( X ) )
		Vepsi_list <- vector( "list", ngroups )
			
		for ( ng in seq( ngroups ) ) {
				
			#- get indices:
			idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
			idx2 <- groupinfo[ng,5]

			#- get data:
			yi    <- y[idx1:idx2]
			Xi    <- X[idx1:idx2,,drop=FALSE] 
			Zpi   <- Zp[idx1:idx2,]
			Zdi   <- Zd[idx1:idx2,]
			Zti   <- Zt[idx1:idx2, ]
			peffi <- person[[ng]]
			deffi <- dyads[[ng]]

			#- compute group-specific Vepsi:
			Vepsi <- if ( !args_list$large ) Zti%*%tcrossprod( G_T, Zti) else tsrm_bayes_compute_Vepsi( L_T=L_T, Zti=Zti, nt=nt)
			
			# note that instead of the usual way: 

			# invVepsi <- base::solve( Vepsi ) 
			# Vepsi_list[[ng]] <- invVepsi
			# tmpXi  <- crossprod( Xi, invVepsi ) # t(Xi)%*%base::solve( Vepsi )
			# sum_mu <- sum_mu + tmpXi%*%( yi - Zpi%*%peffi - Zdi%*%deffi )
			# sum_om <- sum_om + tmpXi%*%Xi
			
			# we use the Cholesky decomposition and backsolve (see technical notes):
			
			Vepsi_chol <- base::chol(Vepsi)
			Vepsi_list[[ng]] <- Vepsi_chol 

			resid <- yi - Zpi %*% peffi - Zdi %*% deffi
			W     <- base::backsolve(Vepsi_chol, cbind(Xi, resid), transpose=TRUE)
			WX    <- W[, seq_len(ncol(Xi)), drop=FALSE]
			Wr    <- W[, ncol(W), drop=FALSE]
			
			#- now the sums:
			sum_mu <- sum_mu + crossprod(WX, Wr)   # = t(Xi)%*%invVepsi%*%resid
			sum_om <- sum_om + crossprod(WX)	   # = t(Xi)%*%invVepsi%*%Xi	

		}

		#- now, get the update:
		OMEGAb <- base::solve( sum_om + invMb )
		MUb    <- OMEGAb %*% ( sum_mu + invMb%*%mb )
		BETA   <- t(mvtnorm::rmvnorm(1,MUb,OMEGAb) )

		#- save all things in a list:
		out <- list( BETA = BETA, OBJ_invVepsi = Vepsi_list )

	}

	#- output:
	return( out )

}