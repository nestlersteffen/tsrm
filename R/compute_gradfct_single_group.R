
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

    #- inverse of V:
	iV <- base::solve( V )
	
	#- residuals + ll-value for ML - part:
	ey    <- as.vector( y - X%*%BETA ) 
	ei    <- crossprod( iV, ey )
	llfct <- mvtnorm::dmvnorm( ey, rep(0,ncol(V)), as.matrix(V), log = TRUE )

	#- derivative for beta:
	dBETA <- crossprod( ei, X )

	#- we pre-compute some things....
	if ( ncol(Zg) > 0 ) { iVZg <- iV %*% Zg; Mg <- crossprod(Zg, iVZg); wg <- crossprod(Zg, ei) }
	iVZp <- iV %*% Zp;  Mp <- crossprod(Zp, iVZp);  wp <- crossprod(Zp, ei)
	iVZd <- iV %*% Zd;  Md <- crossprod(Zd, iVZd);  wd <- crossprod(Zd, ei)
	if (ncol(Zt) > 0) { iVZt <- iV %*% Zt; Mt <- crossprod(Zt, iVZt); wt <- crossprod(Zt, ei) }

	#- REML:
	penvalue_ll <- 0L
	XiVZg <- XiVZp <- XiVZd <- XiVZt <- 0L 
	if ( args_list$with_reml ) {
		XiV <- crossprod(X,iV)
		XiVZg <- if ( ncol(Zg) > 0 ) XiV%*%Zg
		XiVZp <- XiV%*%Zp; XiVZd <- XiV%*%Zd
		XiVZt <- if ( ncol(Zt) > 0 ) XiV%*%Zt
		penvalue_ll <- XiV%*%X 
	}

	#- now compute the gradient for the group:
  	NP   <- max( parm_table$index)
  	NOP  <- nrow( parm_table ) 
  	gradient <- rep( 0, NP )
  	penvalue_gr <- if ( args_list$with_reml ) array( 0, c( ncol(X), ncol(X), NP ) ) else 0L

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
		    	M <- Mp; w <- wp; WZ <- XiVZp
		    } else if ( type %in% c("SD_D","RHO_D") ) { # c("SIGMA_D")
		    	tmp <- diag(1, nd) %x% sigma_derive
		    	tmp <- tmp[perm_d, perm_d]
    			M <- Md; w <- wd; WZ <- XiVZd
		    } else if ( type %in% c("SD_T","RHO_T") ) {  
		    	tmp <- diag(1, nt) %x% sigma_derive
		    	tmp <- tmp[perm_t, perm_t]
			    M <- Mt; w <- wt; WZ <- XiVZt
		    } else if ( type %in% c("SD_G", "RHO_G") ) {
                tmp <- sigma_derive                   
                M <- Mg; w <- wg; WZ <- XiVZg
            }

		    pt1 <- sum(tmp * M)                          # faster version for: tr(P V') = <S, Z^T P Z>
			pt2 <- as.numeric(crossprod(w, tmp %*% w))   # faster version for: tei^T V' tei = w^T S w
		    res <- -0.5*( pt1 - pt2 )
		    
		    #- things for REML penalty term:
		    if ( args_list$with_reml ) {
				penvalue_gr[,,index_nn] <- penvalue_gr[,,index_nn] + WZ%*%tmp%*%t(WZ)
			}

		}

    	gradient[ index_nn ] <- gradient[ index_nn ] + res
    	
  
  	}
	
	return( list( llfct=llfct, gradient=gradient, penvalue_ll=penvalue_ll, penvalue_gr=penvalue_gr ) )
}
