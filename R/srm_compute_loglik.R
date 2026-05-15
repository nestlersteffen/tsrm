
#---- this function computes the log-likelihood for a single group:

srm_compute_loglik_singlegroup <- function( y = NULL, X = NULL, Zg = NULL, Zp = NULL, 
	Zd = NULL, BETA = NULL, SIGMA_G = NULL, SIGMA_P = NULL, SIGMA_D = NULL, 
	random_group = NULL, with_reml = FALSE, fast = FALSE )
{
	#- compute residual:
	ey  <- as.vector( y - X%*%BETA )
	#- compute V:
	V  <- srm_compute_V( Zg, Zp, Zd, SIGMA_G, SIGMA_P, SIGMA_D, random_group ) 
	#- compute the normal distribution:
	res <- mvtnorm::dmvnorm( ey, rep(0,ncol(V)), as.matrix(V), log = TRUE )
	#- add reml?
	if ( with_reml ) {
		iV  <- base::solve( V )
		tmp <- t(X)%*%iV%*%X
		logdet_tmp <- determinant( tmp, logarithm=TRUE )$modulus[1]
		res <- res - 0.5*logdet_tmp
	}
	return( res )
}

#---- this function computes the log-likelihood across groups: 

srm_compute_loglik <- function( parm = NULL, parm_list = NULL, parm_table = NULL,
	data_list = NULL, args_list = NULL, both = NULL ) 
{
	
	#- update parm_list:
  	parm_list <- srm_include_free_parameters( parm = parm, parm_list = parm_list, 
  		parm_table = parm_table )

  	#- get matrices and groupinfo:
  	groupinfo <- data_list[["groupinfo"]]
  	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zg <- data_list[["Zg"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	np <- data_list[["np"]]
	nd <- data_list[["nd"]]
	nv <- data_list[["nv"]]

	#- compute covariance matrices:
	parm_list <- srm_get_sigmas( parm_list = parm_list )

	#- get entries in parm_table and make them full matrices:
	BETA 	<- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]] 
	SIGMA_P <- diag(1,np ) %x% parm_list[["SIGMA_P"]]# %x% diag(1,np ) #
	SIGMA_D <- diag(1,nd ) %x% parm_list[["SIGMA_D"]]# %x% diag(1,nd )  

	#- get no. groups:
	ngroups <- nrow( groupinfo )
	ll <- rep( 0, ngroups )

	#- let's go:
	for ( ng in seq( ngroups ) ) {
		
		#- get indices:
		idx1 <- groupinfo[ng,4] - groupinfo[ng,3] + 1 
		idx2 <- groupinfo[ng,4]
		
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
		
		#- compute loglik - value
		ll[ng] <- srm_compute_loglik_singlegroup( y = y[idx1:idx2], X = X[idx1:idx2,], 
			Zg = Zg[idx1:idx2,], Zp = Zp[idx1:idx2,], Zd = Zd[idx1:idx2,], 
			BETA = BETA_ng, SIGMA_G = SIGMA_G, SIGMA_P = SIGMA_P, SIGMA_D = SIGMA_D, 
			random_group = args_list$random_group, with_reml = args_list$with_reml, 
			fast = args_list$fast )
	}
	
	#- result:
	return( -1*sum(ll) )
}