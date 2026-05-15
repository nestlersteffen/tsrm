
srm_get_ses_beta_singlegroup <- function( parm_list=NULL, parm_table_covs=NULL, 
	np=NULL, nd=NULL, X=NULL, Zg=NULL, Zp=NULL, Zd=NULL, SIGMA_G=NULL, SIGMA_P=NULL, 
	SIGMA_D=NULL, args_list=NULL )
{
	#- compute V and iV:
	V  <- srm_compute_V( Zg=Zg, Zp=Zp, Zd=Zd, SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, 
		SIGMA_D=SIGMA_D, random_group=args_list$random_group )
	iV <- base::solve( V )

	#- compute A and B:
	B <- iV%*%X 
	A <- crossprod( X, B )

	#- Satterthwaite?
	if ( args_list$type_ses == "Satterthwaite") {
		
		#- compute gradient components
  		NP  <- max( parm_table_covs$index)
  		NOP <- nrow( parm_table_covs ) 
  		out <- array( 0, dim = c( ncol( A ), nrow(A), NP + 1 ) )
  		out[,,1] <- A
  	
  	  	for ( nn in 1:NOP ) {
  
		  	#- get matrix type:
	    	free_nn   <- parm_table_covs[nn,]
	    	type      <- free_nn$type
	    	pos       <- c( free_nn$pos1, free_nn$pos2 )
	    	index_nn  <- free_nn$index

	    	#- get derivative:
	    	sigma_derive <- srm_sigma_derivatives( parm_list=parm_list,
	      		type=type, pos=pos )
	    	#- compute derivative of V:
		    if ( type %in% c("SD_P","RHO_P") ) {
		    	tmp <- diag(1,np ) %x% sigma_derive # %x% diag(1,np ) #  
		    	Z   <- Zp
		    } else if ( type %in% c("SD_D","RHO_D") ) { 
		    	tmp <- diag(1,nd ) %x% sigma_derive # %x% diag(1, nd ) #  
		    	Z   <- Zd
		    } else if ( type %in% c("SD_G","RHO_G") ) {  
		    	tmp <- sigma_derive
		    	Z   <- Zg
		    }
		    V_DERIVE <- Z%*%tmp%*%t(Z)
		    out[,,nn+1] <- crossprod( B, V_DERIVE )%*%B
		}

	} else if (args_list$type_ses == "Standard") { 
		out <- array( 0, dim = c( ncol( A ), nrow(A), 1 ) )
		out[,,1] <- A
	}
	#- ...
	return( out )
}

#- function to compute standard errors for the variables:

srm_get_ses <- function( parm=NULL, parm_table=NULL, parm_list=NULL, 
	data_list=NULL, args_list=NULL )  
{

	#- include parm into parm_list:
	parm_list <- srm_include_free_parameters( parm=parm, parm_list=parm_list, 
	  	parm_table=parm_table )

	#- divide parm_table and parm:
	idx <- which( parm_table$type == "BETA")
	parm_table_beta <- parm_table[idx,]
	parm_table_covs <- parm_table[-idx,]
	parm_table_covs$index <- parm_table_covs$index - max( parm_table_beta$index )

	#- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#-  compute standard errors for variance and covariance parms
	tmp_args_list <- args_list
	tmp_args_list$fixed_group <- FALSE
	hessian_covs <- nloptr::nl.jacobian( x0=parm[-idx], fn=srm_compute_gradient, 
	    parm_list=parm_list, parm_table=parm_table_covs, data_list=data_list, 
		args_list=tmp_args_list, both=FALSE )
	vcov_covs <- base::solve( hessian_covs ) 
	ses_covs  <- suppressWarnings( sqrt( diag( vcov_covs ) ) )

	parm_table_covs$est  <- parm[-idx]
	parm_table_covs$se   <- ses_covs
	parm_table_covs$tval <- with( parm_table_covs, est/se )   
	parm_table_covs$p    <- 2 * pnorm( -abs( parm_table_covs$tval ) )

	#- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#-  now we compute the standard errors for the fixed effects
  	groupinfo <- data_list[["groupinfo"]]
  	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zg <- data_list[["Zg"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	np <- data_list[["np"]]
	nd <- data_list[["nd"]]

	#- compute covariance matrices:
	parm_list <- srm_get_sigmas( parm_list = parm_list )

	#- get entries in parm_table and make them full matrices:
	BETA    <- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]]
	SIGMA_P <- parm_list[["SIGMA_P"]] %x% diag(1,np )
	SIGMA_D <- parm_list[["SIGMA_D"]] %x% diag(1,nd )

	#- get no. groups:
	ngroups <- nrow( groupinfo )
	result  <- vector( "list", ngroups )

	#- compute group-specific components:
	for ( ng in seq( ngroups ) ) {

		#- get indices:
		idx1 <- groupinfo[ng,4] - groupinfo[ng,3] + 1 
		idx2 <- groupinfo[ng,4]

		#- compute betas in groups
		result[[ ng ]] <- srm_get_ses_beta_singlegroup( parm_list=parm_list,
			parm_table_covs=parm_table_covs,np=np, nd=nd, X=X[idx1:idx2,], 
			Zg=Zg[idx1:idx2,], Zp=Zp[idx1:idx2,], Zd=Zd[idx1:idx2,], 
			SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, 
			args_list=args_list )
	}

	#- compute sum of matrices:
	result <- Reduce( "+", lapply( result, function( x ) x ) )
	
	#- compute A and its inverse
	A  <- result[,,1]
	vcov_beta <- base::solve( A )
	ses_beta  <- sqrt( diag( vcov_beta ) )
	parm_table_beta$est  <- parm[idx]
	parm_table_beta$se   <- ses_beta
	parm_table_beta$tval <- with( parm_table_beta, est/se )
	
	#- for the p-value we use the standard test or the Sattertwhaite approximation
	if ( args_list$type_ses == "Satterthwaite" ) {
		p    <- ncol(X)
		cmat <- diag(1,p) 
	   	dfs  <- rep(0,p)
	   	for (pp in seq(p) ) {
	      	cvec <- cmat[,pp]
	      	#- get gradient
	      	gr <- rep( 0, dim( result )[3] - 1 )
			for ( gg in seq( length( gr ) ) ) {
	      		gr[gg] <- (t(cvec)%*%(vcov_beta%*%result[,,gg+1]%*%vcov_beta))%*%cvec 
	      	}
	      	#- get approximate dfs
	      	dfs[pp] <- (2*(t(cvec)%*%vcov_beta%*%cvec)^2)/(t(gr)%*%vcov_covs%*%gr)
	   } 
	   parm_table_beta$dfs <- dfs
	   parm_table_beta$p   <- 2*pt(-abs( parm_table_beta$tval), df=abs(dfs) )
	   parm_table_covs$dfs <- rep( 1, nrow( parm_table_covs) )
	} else if ( args_list$type_ses == "Standard" ) { 
		parm_table_beta$p  <- 2 * pnorm( -abs( parm_table_beta$tval ) )  
	}

	#- combine tables and add final information:
	parm_table_covs$index <- parm_table_covs$index + max( parm_table_beta$index )
	parm_table <- rbind( parm_table_beta, parm_table_covs )
	parm_table$lCI <- parm_table$est - 1.96*parm_table$se
	parm_table$uCI <- parm_table$est + 1.96*parm_table$se 
	
	#- return parm_table
	return( parm_table )
}