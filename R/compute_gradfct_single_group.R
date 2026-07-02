
#---- this function computes the log-likelihood and the gradient for a single group:

compute_gradient_singlegroup <- function( parm_list=NULL, parm_table=NULL, 
	np=NULL, nd=NULL, nt=NULL, 
	y=NULL, X=NULL, Zg=NULL, Zp=NULL, Zd=NULL, Zt=NULL, 
	perm_p=NULL, perm_d=NULL, perm_t=NULL,
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

	#- we pre-compute some things....
	if (ncol(Zg) > 0) { PZg <- P %*% Zg; Mg <- crossprod(Zg, PZg); wg <- crossprod(Zg, ei) }
	PZp <- P %*% Zp;  Mp <- crossprod(Zp, PZp);  wp <- crossprod(Zp, ei)
	PZd <- P %*% Zd;  Md <- crossprod(Zd, PZd);  wd <- crossprod(Zd, ei)
	if (ncol(Zt) > 0) { PZt <- P %*% Zt; Mt <- crossprod(Zt, PZt); wt <- crossprod(Zt, ei) }

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
		    	tmp <- diag(1, np) %x% sigma_derive
		    	tmp <- tmp[perm_p, perm_p] # faster than: tmp <- Pp %*% ( diag(1, np ) %x% sigma_derive ) %*% t( Pp )
		    	M <- Mp; w <- wp
		    } else if ( type %in% c("SD_D","RHO_D") ) { # c("SIGMA_D")
		    	tmp <- diag(1, nd) %x% sigma_derive
		    	tmp <- tmp[perm_d, perm_d]
    			M <- Md; w <- wd
		    } else if ( type %in% c("SD_T","RHO_T") ) {  
		    	tmp <- diag(1, nt) %x% sigma_derive
		    	tmp <- tmp[perm_t, perm_t]
			    M <- Mt; w <- wt
		    } else if ( type %in% c("SD_G", "RHO_G") ) {
                tmp <- sigma_derive                   
                M <- Mg; w <- wg
            }

		    pt1 <- sum(tmp * M)                          # faster version for: tr(P V') = <S, Z^T P Z>
			pt2 <- as.numeric(crossprod(w, tmp %*% w))   # faster version for: tei^T V' tei = w^T S w
		    res 	 <- -0.5*( pt1 - pt2 )

		}

    	ng_grad[ index_nn ] <- ng_grad[ index_nn ] + res
  
  	}
	
	return( c( ng_ll, ng_grad ) )
}
