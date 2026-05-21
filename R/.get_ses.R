
get_ses_beta_singlegroup <- function( parm_list=NULL, parm_table_covs=NULL, 
	np=NULL, nd=NULL, nt=NULL, 
	X=NULL, Zg=NULL, Zp=NULL, Zd=NULL, Zt=NULL, 
	SIGMA_G=NULL, SIGMA_P=NULL, SIGMA_D=NULL, SIGMA_T=NULL, 
	args_list=FALSE, sigma_derivatives=NULL )
{
	#- compute V and iV:
    V <- Zp %*% SIGMA_P %*% t(Zp) +
         Zd %*% SIGMA_D %*% t(Zd)

    if ( ncol(Zt) > 0 ) {
        V <- V + Zt %*% SIGMA_T %*% t(Zt)
    }
    if ( ncol(Zg) > 0 ) {
        V <- V + Zg %*% SIGMA_G %*% t(Zg)
    }
	iV <- base::solve( V )

	#- compute A and B:
	B <- iV%*%X 
	A <- crossprod( X, B )

	#- Satterthwaite?
	if ( args_list$type_ses == "Satterthwaite") {
		
		#- compute gradient components
  		NP  	 <- max( parm_table_covs$index)
  		NOP 	 <- nrow( parm_table_covs ) 
  		out      <- array( 0, dim = c( ncol( A ), nrow(A), NP + 1 ) )
  		out[,,1] <- A
  	
  	  	for ( nn in 1:NOP ) {
  
		  	#- get matrix type:
	    	free_nn   <- parm_table_covs[nn,]
	    	type      <- free_nn$type
	    	pos       <- c( free_nn$pos1, free_nn$pos2 )
	    	index_nn  <- free_nn$index

	    	#- get derivative:
	 		sigma_derive <- sigma_derivatives( parm_list=parm_list,
	      		type=type, pos=pos )
	  
		    #- compute derivative of V depending on matrix:
		    if ( type %in% c("SD_P","RHO_P") ) { # c("SIGMA_P")
		    	tmp <- diag(1, np ) %x% sigma_derive
		    	Z   <- Zp
		    } else if ( type %in% c("SD_D","RHO_D") ) { # c("SIGMA_D")
		    	tmp <- diag(1, nd ) %x% sigma_derive
		    	Z   <- Zd
		    } else if ( type %in% c("SD_T","RHO_T") ) {  
		    	tmp <- diag(1, nt ) %x% sigma_derive
		    	Z   <- Zt
		    } else if ( type %in% c("SD_G", "RHO_G") ) {
                tmp <- sigma_derive                   
                Z   <- Zg
            }

		    V_DERIVE <- Z%*%tmp%*%t(Z)
		    out[,,nn+1] <- crossprod( B, V_DERIVE )%*%B
		}

	} else if ( args_list$type_ses == "Standard" ) { 
		out <- array( 0, dim = c( ncol( A ), nrow(A), 1 ) )
		out[,,1] <- A
	}
	
	return( out )
}

#- function to compute standard errors for the variables:

get_ses <- function( parm=NULL, parm_table=NULL, parm_list=NULL, data_list=NULL, 
	args_list=NULL, model=c("srm","tsrm") )  
{

	#- get model-specific functions.
	model <- match.arg( model )
	if ( model == "srm" ) {
		include_free_parms <- srm_include_free_parameters
        sigma_derivatives  <- srm_sigma_derivatives
	} else {
		include_free_parms <- tsrm_include_free_parameters
        sigma_derivatives  <- tsrm_sigma_derivatives
	}

	#- insert parm
	parm_list <- include_free_parms( parm=parm, parm_list=parm_list, 
		parm_table=parm_table )

	#- divide parm_table and parm:
	idx <- which( parm_table$type == "BETA")
	parm_table_beta <- parm_table[idx,]
	parm_table_covs <- parm_table[-idx,]
	parm_table_covs$index <- parm_table_covs$index - max( parm_table_beta$index )

	#- compute standard errors for variance and covariance parms:
	hessian_covs <- nloptr::nl.jacobian( x0=parm[-idx], fn=compute_gradfct, 
	    parm_list=parm_list, parm_table=parm_table_covs, data_list=data_list, 
		args_list=args_list, model=model, both=FALSE )
	vcov_covs 			<- base::solve( hessian_covs ) 
	ses_covs  			<- suppressWarnings( sqrt( diag( vcov_covs ) ) )
	parm_table_covs$est <- parm[-idx]
	parm_table_covs$se  <- ses_covs
	parm_table_covs$z   <- with( parm_table_covs, est/se )   
	parm_table_covs$p   <- 2 * pnorm( -abs( parm_table_covs$z ) )

	#- now we compute the standard errors for the fixed effects:
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
	
	#- compute covariance matrices:
	parm_list <- get_sigmas( parm_list=parm_list, model=model  )

	#- get entries in parm_table and make them full matrices:
	BETA    <- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]]
	SIGMA_P <- parm_list[["SIGMA_P"]] %x% diag(1,np )
	SIGMA_D <- parm_list[["SIGMA_D"]] %x% diag(1,nd )
	SIGMA_T <- if ( ncol(Zt) > 0 ) diag(1, nt) %x% parm_list[["SIGMA_T"]] else parm_list[["SIGMA_T"]]

	#- get no. groups:
	ngroups <- nrow( groupinfo )
	result  <- vector( "list", ngroups )

	#- compute group-specific components:
	for ( ng in seq( ngroups ) ) {

		#- get indices:
		idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
		idx2 <- groupinfo[ng,5]

		#- get group and triadic matrices:
		Zg_ng <- if (ncol(Zg) > 0) as.matrix( Zg[idx1:idx2, ] ) else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)
		Zt_ng <- if (ncol(Zt) > 0) Zt[idx1:idx2, ] else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)

		#- compute betas in groups
		result[[ ng ]] <- get_ses_beta_singlegroup( parm_list=parm_list,
			parm_table_covs=parm_table_covs,np=np, nd=nd, nt=nt, X=X[idx1:idx2,], 
			Zg=Zg_ng, Zp=Zp[idx1:idx2,], Zd=Zd[idx1:idx2,], Zt=Zt_ng,
			SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, SIGMA_T=SIGMA_T, 
			args_list=args_list, sigma_derivatives=sigma_derivatives )
	}

	#- compute sum of matrices:
	result <- Reduce( "+", lapply( result, function( x ) x ) )
	
	#- compute A and its inverse:
	A <- result[,,1]
	vcov_beta <- base::solve( A )
	ses_beta  <- sqrt( diag( vcov_beta ) )
	parm_table_beta$est <- parm[idx]
	parm_table_beta$se  <- ses_beta
	parm_table_beta$z   <- with( parm_table_beta, est/se )
	
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
	   parm_table_beta$p   <- 2*pt(-abs( parm_table_beta$z), df=abs(dfs) )
	   parm_table_covs$dfs <- rep( 1, nrow( parm_table_covs) )
	
	} else if ( args_list$type_ses == "Standard" ) { 
	
		parm_table_beta$p  <- 2 * pnorm( -abs( parm_table_beta$z ) )  
	
	}

	#- combine tables and add final information:
	parm_table_covs$index <- parm_table_covs$index + max( parm_table_beta$index )
	parm_table <- rbind( parm_table_beta, parm_table_covs )
	parm_table$lCI <- parm_table$est - 1.96*parm_table$se
	parm_table$uCI <- parm_table$est + 1.96*parm_table$se 
	
	#- return parm_table
	return( parm_table )
}