
#---- this function computes the log-likelihood and the gradient for a single group:

compute_gradient_singlegroup <- function( parm_list=NULL, parm_table=NULL, 
	np=NULL, nd=NULL, nt=NULL, 
	y=NULL, X=NULL, Zg=NULL, Zp=NULL, Zd=NULL, Zt=NULL, Pp=NULL, Pd=NULL, 
	BETA=NULL, SIGMA_G=NULL, SIGMA_P=NULL, SIGMA_D=NULL, SIGMA_T=NULL, 
	args_list=FALSE, sigma_derivatives=NULL )
{
	
	#- compute V:
    V <- Zp %*% SIGMA_P %*% t(Zp) +
         Zd %*% SIGMA_D %*% t(Zd)

    if ( ncol(Zt) > 0 ) {
        V <- V + Zt %*% SIGMA_T %*% t(Zt)
    }
    if ( ncol(Zg) > 0 ) {
        V <- V + Zg %*% SIGMA_G %*% t(Zg)
    }

	iV <- base::solve( V )
	P  <- iV
	
	#- residuals + ll:
	ey    <- as.vector( y - X%*%BETA ) 
	ei    <- crossprod( iV, ey )
	ng_ll <- mvtnorm::dmvnorm( ey, rep(0,ncol(V)), as.matrix(V), log = TRUE )

	#- REML:
	if ( args_list$with_reml ) {
		tmp   <- crossprod(X,iV)%*%X
		ng_ll <- ng_ll - 0.5*determinant( tmp, logarithm=TRUE )$modulus[1]
		P     <- iV - (iV%*%X)%*%base::solve(tmp)%*%(crossprod(X,iV))
	}
	
	#- derivative for beta:
	dBETA <- crossprod( ei, X )
	
	#- now compute the gradient for the group:
  	NP      <- max( parm_table$index)
  	NOP     <- nrow( parm_table ) 
  	ng_grad <- rep( 0, NP )

  	for ( nn in 1:NOP ) {
  
	  	#- get matrix type:
	    free_nn   <- parm_table[nn,]
	    type      <- free_nn$type
	    pos       <- c( free_nn$pos1, free_nn$pos2 )
	    index_nn  <- free_nn$index

	    if ( type == "BETA" ) {

			res <- dBETA[ pos[1] ] 		
 		
 		} else {
 	
	 		#- get derivative:
	 		sigma_derive <- sigma_derivatives( parm_list=parm_list,
	      		type=type, pos=pos )
	  
		    #- compute derivative of V depending on matrix:
		    if ( type %in% c("SD_P","RHO_P") ) { # c("SIGMA_P")
		    	tmp <- Pp %*% ( diag(1, np ) %x% sigma_derive ) %*% t( Pp )
		    	Z   <- Zp
		    } else if ( type %in% c("SD_D","RHO_D") ) { # c("SIGMA_D")
		    	tmp <- Pd %*% ( diag(1, nd ) %x% sigma_derive ) %*% t( Pd )
		    	Z   <- Zd
		    } else if ( type %in% c("SD_T","RHO_T") ) {  
		    	tmp <- diag(1, nt ) %x% sigma_derive
		    	Z   <- Zt
		    } else if ( type %in% c("SD_G", "RHO_G") ) {
                tmp <- sigma_derive                   
                Z   <- Zg
            }

		    V_DERIVE <- Z%*%tmp%*%t(Z)
		    pt1      <- sum( P * V_DERIVE )
		    pt2 	 <- ( t( ei ) %*% V_DERIVE ) %*% ei 
		    res 	 <- -0.5*( pt1 - pt2 )

		}

    	ng_grad[ index_nn ] <- ng_grad[ index_nn ] + res
  
  	}
	
	return( c( ng_ll, ng_grad ) )
}
