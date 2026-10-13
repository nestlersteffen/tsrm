
#---- function to update random person and dyad effects
#     Note: function uses pointer or list of inverses of Veps, computed in BETA update

tsrm_bayes_update_effects <- function( data_list=NULL, args_list=NULL, person=NULL, dyads=NULL, 
	BETA=NULL, SIGMA_P=NULL, SIGMA_D=NULL, OBJ_invVepsi=NULL ) 
{

	#- get matrices and groupinfo:
  	groupinfo <- data_list[["groupinfo"]]
  	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	np <- data_list[["np"]]
	nd <- data_list[["nd"]]
	nt <- data_list[["nt"]]

	#- get no. groups:
	ngroups <- nrow( groupinfo )

	#- compute relevant matrices:
	invG_P <- diag( 1, np ) %x% base::solve( SIGMA_P )
	invG_D <- diag( 1, nd ) %x% base::solve( SIGMA_D )
	
	if ( args_list$use_rcpp ) {

		#- change groupinfo to integer
		storage.mode( groupinfo ) <- "integer"

		#- call the C++ function
		out <- tsrm_bayes_update_effects_export( y = y, X = X, Zp = Zp, Zd = Zd,  
			invG_P = invG_P, invG_D = invG_D, BETA = as.vector( BETA ), 
			groupinfo = groupinfo, ngroups = ngroups, person = person, 
			dyads = dyads, OBJ_invVepsi = OBJ_invVepsi )
	
	} else {

		#- make a list with person effects
		new_person <- vector( "list", ngroups )
		new_dyads  <- vector( "list", ngroups )

		for ( ng in seq( ngroups ) ) {
				
			#- get indices:
			idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
			idx2 <- groupinfo[ng,5]

			#- get data:
			yi    <- y[idx1:idx2]
			Xi    <- X[idx1:idx2,] 
			Zpi   <- Zp[idx1:idx2,]
			Zdi   <- Zd[idx1:idx2,]
			peffi <- person[[ng]]
			deffi <- dyads[[ng]]

			#- compute group-specific Veps:
			invVepsi <- OBJ_invVepsi[[ ng ]]

			#- update person effects OMEGAp <- (OMEGAp + t(OMEGAp)) / 2  # Symmetrie erzwingen
			tmpZVp <- crossprod( Zpi, invVepsi )
			OMEGAp <- base::solve( tmpZVp %*% Zpi + invG_P )
			MUp	   <- OMEGAp %*% ( tmpZVp %*%( yi - Xi%*%BETA - Zdi%*%deffi ) )
			new_peffi <- t( mvtnorm::rmvnorm( 1, MUp, as.matrix( OMEGAp ) ) )
				
			#- update dyad effects (with new person-effects):
			tmpZVd <- crossprod( Zdi, invVepsi )
			OMEGAd <- base::solve( tmpZVd %*% Zdi + invG_D )
			MUd	   <- OMEGAd %*% ( tmpZVd %*%( yi - Xi%*%BETA - Zpi%*%new_peffi ) )
			new_deffi <- t( mvtnorm::rmvnorm( 1, MUd, as.matrix( OMEGAd ) ) )
				
			#- save in list:
			new_person[[ng]] <- new_peffi
			new_dyads[[ng]]  <- new_deffi

		}

		#- output:
		out <- list( person = new_person, dyads = new_dyads )
	
	}

	return( out )

}