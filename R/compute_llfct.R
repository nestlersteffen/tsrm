
#---- this function computes the log-likelihood for all groups:

compute_llfct <- function( parm=NULL, parm_list=NULL, parm_table=NULL,
	data_list=NULL, args_list=NULL, model=c("srm","tsrm") )
{

	#- get model-specific functions.
	model <- match.arg( model )
	if ( model == "srm" ) {
		include_free_parms <- srm_include_free_parameters
	} else {
		include_free_parms <- tsrm_include_free_parameters
	}

	#- insert parm
	parm_list <- include_free_parms( parm=parm, parm_list=parm_list, 
		parm_table=parm_table )

	#- get covariance matrices
	parm_list <- get_sigmas( parm_list=parm_list, model=model )

	#- get parm matrices:
	BETA 	<- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]]
	SIGMA_P <- parm_list[["SIGMA_P"]]
	SIGMA_D <- parm_list[["SIGMA_D"]]
    SIGMA_T <- parm_list[["SIGMA_T"]]

	#- get matrices and groupinfo:
	groupinfo <- data_list[["groupinfo"]]
	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	Zg <- data_list[["Zg"]] 
	Zt <- data_list[["Zt"]] 
	np <- data_list[["np"]]
	nd <- data_list[["nd"]]
	nt <- data_list[["nt"]]
	nv <- data_list[["nv"]]
	Pp <- data_list[["Pp"]]
	Pd <- data_list[["Pd"]]
	
	#- make the matrices "big"
	SIGMA_P <- Pp %*% ( diag(1,np) %x% SIGMA_P ) %*% t( Pp )
	SIGMA_D <- Pd %*% ( diag(1,nd) %x% SIGMA_D ) %*% t( Pd )
	SIGMA_T <- if ( ncol(Zt) > 0 ) diag(1, nt) %x% parm_list[["SIGMA_T"]] else SIGMA_T

	#- get no. groups:
	ngroups  <- nrow( groupinfo )
	llfct    <- rep( 0, ngroups )

	for ( ng in seq( ngroups ) ) {
				
		#- get indices:
		idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
		idx2 <- groupinfo[ng,5]

		#- fixed groups?
		BETA_ng <- BETA
		if ( args_list$fixed_group ) {
			if ( length( BETA ) > ngroups ) {
				idx <- c( ng, (ngroups+1):length(BETA) )
			} else {
				idx <- ng
			}
			BETA_ng <- as.matrix( BETA[idx,1] )
		}

		#- get group specific matrices:
		Zg_ng <- if (ncol(Zg) > 0) as.matrix( Zg[idx1:idx2, ] ) else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)
		Zp_ng <- Zp[idx1:idx2,]
		Zd_ng <- Zd[idx1:idx2,]
		Zt_ng <- if (ncol(Zt) > 0) Zt[idx1:idx2, ] else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)
		X_ng  <- X[idx1:idx2,]
		y_ng  <- y[idx1:idx2] 		

		#- compute mean:
		Xb <- X_ng%*%BETA_ng

		#- compute covariance matrix V:
		V <- Zp %*% SIGMA_P %*% t(Zp) +
			 Zd %*% SIGMA_D %*% t(Zd)

    	if ( ncol(Zt) > 0 ) {
        	V <- V + Zt %*% SIGMA_T %*% t(Zt)
    	}
    	if ( ncol(Zg) > 0 ) {
        	V <- V + Zg %*% SIGMA_G %*% t(Zg)
    	}

		#- compute loglik - value
		res <- mvtnorm::dmvnorm( y_ng, Xb, as.matrix(V), log=TRUE )

		#- add REML:
		if ( args_list$with_reml ) {
			iV  <- base::solve( V )
			tmp <- t(X_ng)%*%iV%*%X_ng
			logdet_tmp <- determinant( tmp, log = TRUE )$modulus[1]
			res <- res - 0.5*logdet_tmp
		}

		llfct[ ng ] <- res

	}
	
	#- result:
	return( -1*sum( llfct ) )

}