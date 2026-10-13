
#---- function to compute the update of BETA in TSRM-Bayes 

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
	G_T <- diag( 1, nt ) %x% SIGMA_T

	#- get hyperparms:
	mb 	  <- args_list$mcmc_args$mb
	invMb <- args_list$mcmc_args$invMb

	if ( args_list$use_rcpp ) {

		#- change groupinfo to integer:
		storage.mode( groupinfo ) <- "integer"

		#- call the C++ function
		tmp <- tsrm_bayes_update_beta_export( y=y, X=X, Zp=Zp, Zd=Zd, Zt=Zt, 
			G_T=G_T, groupinfo=groupinfo, ngroups=ngroups, person=person, 
			dyads=dyads, mb=mb, invMb=invMb )
		out <- list( BETA = as.matrix( tmp$BETA ), OBJ_invVepsi = tmp$OBJ_invVepsi )

	} else {

		#- initialize:
		sum_mu <- matrix( 0, nrow = ncol( X ), ncol = 1 )
		sum_om <- diag( 0, ncol( X ) )
		invVepsi_list <- vector( "list", ngroups )
			
		for ( ng in seq( ngroups ) ) {
				
			#- get indices:
			idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
			idx2 <- groupinfo[ng,5]

			#- get data:
			yi    <- y[idx1:idx2]
			Xi    <- X[idx1:idx2,] 
			Zpi   <- Zp[idx1:idx2,]
			Zdi   <- Zd[idx1:idx2,]
			Zti   <- Zt[idx1:idx2, ]
			peffi <- person[[ng]]
			deffi <- dyads[[ng]]

			#- compute group-specific Veps:
			Vepsi <- Zti%*%tcrossprod( G_T, Zti) # Zdi%*%G_D%*%t(Zdi)
			invVepsi <- base::solve( Vepsi )
			invVepsi_list[[ng]] <- invVepsi

			#- now the sums:
			tmpXi  <- crossprod( Xi, invVepsi ) # t(Xi)%*%base::solve( Vepsi )
			sum_mu <- sum_mu + tmpXi%*%( yi - Zpi%*%peffi - Zdi%*%deffi )
			sum_om <- sum_om + tmpXi%*%Xi

		}

		#- now, get the update:
		OMEGAb <- base::solve( sum_om + invMb )
		MUb    <- OMEGAb %*% ( sum_mu + invMb%*%mb )
		BETA   <- t(mvtnorm::rmvnorm(1,MUb,OMEGAb) )

		#- save all things in a list:
		out <- list( BETA = BETA, OBJ_invVepsi = invVepsi_list )

	}

	#- output:
	return( out )

}

# faster?
# cholVeps <- chol(Vepsi)
        
# # Compute residuals
# resid <- yi - Zpi %*% peffi
        
# # Solve linear systems using Cholesky
# # Instead of solve(Vepsi, resid), use:
# tmp_resid <- backsolve(cholVeps, backsolve(cholVeps, resid, transpose = TRUE))
# tmp_Xi <- backsolve(cholVeps, backsolve(cholVeps, Xi, transpose = TRUE))
        
# # Accumulate sufficient statistics
# sum_mu <- sum_mu + crossprod(Xi, tmp_resid)  # Faster than t(Xi) %*% ...
# sum_om <- sum_om + crossprod(Xi, tmp_Xi)