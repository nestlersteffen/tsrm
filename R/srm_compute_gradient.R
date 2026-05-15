
#---- this function computes the log-likelihood for a single group:

srm_compute_gradient_singlegroup <- function( parm_list=NULL, parm_table=NULL, 
	np=NULL, nd=NULL, y=NULL, X=NULL, Zg=NULL, Zp=NULL, Zd=NULL, 
	BETA=NULL, SIGMA_G=NULL, SIGMA_P=NULL, SIGMA_D=NULL, 
	random_group=NULL, with_reml=FALSE )
{
	
	#- compute V or inverse of V and adapt P:
	V  <- srm_compute_V( Zg, Zp, Zd, SIGMA_G, SIGMA_P, SIGMA_D, random_group=random_group ) 
	iV <- base::solve( V )
	P  <- iV

	#- compute residuals:
	ey <- as.vector( y - X%*%BETA ) 
	ei <- crossprod( iV, ey )

	#- compute log-likelihood:
	ng_ll <- mvtnorm::dmvnorm( ey, rep(0,ncol(V)), as.matrix(V), log = TRUE )

	#- compute REML part:
	if ( with_reml ) {
		tmp   <- crossprod(X,iV)%*%X
		ng_ll <- ng_ll - 0.5*determinant( tmp, logarithm=TRUE )$modulus[1]
		P     <- iV - (iV%*%X)%*%base::solve(tmp)%*%(crossprod(X,iV))
	}
	
	#- derivative for beta:
	dBETA <- crossprod( ei, X )
	
	#- now compute the gradient for the group:
  	NP  <- max( parm_table$index)
  	NOP <- nrow( parm_table ) 
  	ng_grad <- rep( 0, NP )
  	
  	for ( nn in 1:NOP ) {
  
	  	#- get matrix type:
	    free_nn   <- parm_table[nn,]
	    type      <- free_nn$type
	    pos       <- c( free_nn$pos1, free_nn$pos2 )
	    index_nn  <- free_nn$index

	    if ( type == "BETA" ) {
			
			res <- dBETA[ pos[1] ]
 			V_DERIVE <- NULL
 		
 		} else {
 	
	 		#- get derivative:
	 		sigma_derive <- srm_sigma_derivatives( parm_list = parm_list,
	      		type = type, pos = pos )
	  
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
		  	
		    #- compute gradient value
		    pt1 <- sum( P * V_DERIVE )
		    pt2 <- ( t( ei ) %*% V_DERIVE ) %*% ei 
		    res <- -0.5*( pt1 - pt2 )

		}
    
    	ng_grad[ index_nn ] <- ng_grad[ index_nn ] + res
  
  	}
	
	#- output:
	return( c( ng_ll, ng_grad ) )
}

#---- this function computes the log-likelihood across groups: 

srm_compute_gradient <- function( parm = NULL, parm_list = NULL, parm_table = NULL,
	data_list = NULL, args_list = NULL, both = FALSE ) 
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

	#- compute covariance matrices:
	parm_list <- srm_get_sigmas( parm_list = parm_list )

	#- get entries in parm_table and make them full matrices:
	BETA 	<- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]]
	SIGMA_P <- diag(1,np ) %x% parm_list[["SIGMA_P"]] # %x% diag(1,np ) #
	SIGMA_D <- diag(1,nd ) %x% parm_list[["SIGMA_D"]] # %x% diag(1,nd ) # 

	#- get no. groups:
	ngroups  <- nrow( groupinfo )
	gradient <- matrix( 0, nrow = ngroups, ncol = nrow( parm_table$index) ) 
	llfct    <- rep( 0, ngroups )

	#- let's go:
	for ( ng in seq( ngroups ) ) {
		
		#- get indices:
		idx1 <- groupinfo[ng,4] - groupinfo[ng,3] + 1 
		idx2 <- groupinfo[ng,4]
		
		#- fixed groups?
		BETA_ng <- BETA
		parm_table_ng <- parm_table
		if ( args_list$fixed_group ) {
			if ( length( BETA ) > ngroups ) {
				idx <- c( ng, (ngroups+1):length(BETA) )
			} else {
				idx <- ng 
			}
			BETA_ng <- as.matrix( BETA[idx,1] )
			parm_table_ng <- parm_table[ (parm_table$type == "BETA" & parm_table$pos1 %in% idx) | 
				parm_table$type != "BETA", ]
			parm_table_ng[ parm_table_ng$type == "BETA", ]$pos1 <- seq( 1, nrow( BETA_ng ), 1 )
		}

		#- compute all relevant values
		if ( args_list$use_rcpp && !args_list$large ) {

			#- get position matrix and type in case of rcpp:
			posmat  <- as.matrix( parm_table_ng[,c("pos1","pos2")] )
			typevec <- as.vector( parm_table_ng[,c("ntype")] )
			storage.mode( typevec ) <- "integer"
			storage.mode( posmat )  <- "integer" 
			
			#- compute ...
			tmp_grad <- srm_gradient_singlegroup_export( parm_list=parm_list, 
				ty=as.matrix(y[idx1:idx2]), tX=as.matrix(X[idx1:idx2,]), 
				tZg=as.matrix(Zg[idx1:idx2,]), 
				tZp=Zp[idx1:idx2,], tZd=Zd[idx1:idx2,], 
				SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, 
				BETA=as.vector(BETA_ng), np=np, nd=nd, typevec=typevec, posmat=posmat, 
				with_reml=args_list$with_reml, random_group=args_list$random_group )

		} else if ( args_list$use_rcpp && args_list$large ) {

			#- make sparse matrices:
			tZp <- as(Zp[idx1:idx2,], "dgCMatrix")
			tZd <- as(Zd[idx1:idx2,], "dgCMatrix")
			SIGMA_P <- as(SIGMA_P, "dgCMatrix")
			SIGMA_D <- as(SIGMA_D, "dgCMatrix")

			#- get position matrix and type in case of rcpp:
			posmat  <- as.matrix( parm_table_ng[,c("pos1","pos2")] )
			typevec <- as.vector( parm_table_ng[,c("ntype")] )
			storage.mode( typevec ) <- "integer"
			storage.mode( posmat )  <- "integer" 
			
			#- compute ...
			tmp_grad <- srm_gradient_singlegroup_sparse_export( parm_list=parm_list, 
				ty=as.matrix(y[idx1:idx2]), tX=as.matrix(X[idx1:idx2,]), 
				tZg=as.matrix(Zg[idx1:idx2,]), 
				tZp=tZp, tZd=tZd, 
				SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, 
				BETA=as.vector(BETA_ng), np=np, nd=nd, typevec=typevec, posmat=posmat, 
				with_reml=args_list$with_reml, random_group=args_list$random_group )

		} else {
			#- compute ...
			tmp_grad <- srm_compute_gradient_singlegroup( parm_list=parm_list, 
				parm_table=parm_table_ng, np=np, nd=nd, y=y[idx1:idx2], X=X[idx1:idx2,], 
				Zg=Zg[idx1:idx2,], Zp=Zp[idx1:idx2,], Zd=Zd[idx1:idx2,], BETA=BETA_ng, 
				SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, 
				random_group=args_list$random_group, with_reml=args_list$with_reml )
		}
		
		llfct[ng]     <- tmp_grad[ 1]
		gradient[ng,] <- tmp_grad[-1]
		
	}
	#- output:
	if ( both ) {
		return( list( "objective"=-1*sum(llfct), "gradient"=-1*colSums(gradient) ) )
	} else {
		return( -1*colSums(gradient) )
	}
}