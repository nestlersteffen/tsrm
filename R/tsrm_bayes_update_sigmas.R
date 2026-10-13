

tsrm_bayes_update_sigmap <- function( person=NULL, args_list=NULL, sum_np=NULL, old=NULL ) 
{

	#- compute residual sum of square matrix ( sum p*p^T )
	MAT <- matrix( do.call( "c", person ), ncol=3, byrow=TRUE ) 
	RES <- t( MAT )%*%( MAT )

	#- which prior?
	prior_sigmap <- args_list$mcmc_args$prior_sigmap

	if ( prior_sigmap == "IW" ) {

		#- get hyperparms:
		nup <- args_list$mcmc_args$IW.nup 
		Sp  <- args_list$mcmc_args$IW.Sp 
			
		#- generate draw:
		SIGMA_P <- MCMCpack::riwish( nup + sum_np, Sp + RES )
		R_P     <- cov2cor( SIGMA_P )
		update  <- c( sqrt( diag( SIGMA_P ) ), R_P[lower.tri( R_P, diag = FALSE)] )
		accept  <- 1
	
	} else if ( prior_sigmap == "LKJ" ) {
		
		#- collect hyperparms:
		eta   <- args_list$mcmc_args$LKJ.eta
		scale <- args_list$mcmc_args$Ht.scale
		nu    <- args_list$mcmc_args$Ht.nu
		tau   <- args_list$mcmc_args$tau

		#- make the matrices for the current proposal:
		c_sd      <- c( old[1], old[2],old[3] )
		c_SD_P    <- diag( c_sd, 3 )
		c_RHO_P   <- matrix(c( 1,old[4],old[5],old[4],1,old[6],old[5],old[6],1),3,3)
		c_SIGMA_P <- c_SD_P %*% c_RHO_P %*% c_SD_P

		#- generate new standard deviations and new correlation matrix:
		p_logsd <- rnorm( 3, mean=log( c( old[1], old[2],old[3] ) ), sd=tau )
		p_sd    <- exp( p_logsd )
		p_SD_P  <- diag( p_sd, 3 )
		p_RHO_P <- bayes_rlkjcorr( 1, 3, eta=eta )

		#- generate covariance matrix:
		p_SIGMA_P <- p_SD_P%*%p_RHO_P%*%p_SD_P

		#- compute acceptance probability
		cLogPosterior  <- -(sum_np/2)*base::determinant( c_SIGMA_P, log = TRUE )$modulus[1] -
			0.5*sum( RES * base::solve( c_SIGMA_P ) ) +
			sum( bayes_dhalft( c_sd, scale=scale, nu=nu, log=TRUE ) ) +
			bayes_dlkjcorr( c_RHO_P, eta=eta, log=TRUE )
		pLogPosterior  <- -(sum_np/2)*base::determinant( p_SIGMA_P, log = TRUE )$modulus[1] -
			0.5*sum( RES * base::solve( p_SIGMA_P ) ) +
			sum( bayes_dhalft( p_sd, scale=scale, nu=nu, log=TRUE ) ) +
			bayes_dlkjcorr( p_RHO_P, eta=eta, log=TRUE )
		cLogProposal   <- bayes_dlkjcorr( c_RHO_P, eta=eta, log=TRUE )
		pLogProposal   <- bayes_dlkjcorr( p_RHO_P, eta=eta, log=TRUE ) 
		LogJacobiRatio <- sum( p_logsd ) - sum( log( c_sd ) )	
		llratio <- (pLogPosterior - cLogPosterior) + ( cLogProposal - pLogProposal ) + LogJacobiRatio

		#- now accept or reject both together
		if ( log(runif(1)) < llratio ) {
	    	SIGMA_P <- p_SIGMA_P
	    	update <- c( p_sd, p_RHO_P[ rbind( c(1,2),c(1,3),c(2,3) ) ] )
	    	accept <- 1  
		} else {
			SIGMA_P <- c_SIGMA_P
			update <- old
			accept <- 0  
		}

	} else if ( prior_sigmap == "LKJ_single" ) {

		#- hyperparameters:
		eta       <- args_list$mcmc_args$LKJ.eta
		scale     <- args_list$mcmc_args$Ht.scale
		nu        <- args_list$mcmc_args$Ht.nu
		tau_sd    <- args_list$mcmc_args$tau_sd    # RW-sd for log(sd), e.g. 0.1-0.2
		tau_rho   <- args_list$mcmc_args$tau_rho   # RW-sd for logit(u), e.g. 0.3-0.6

		#- current state:
		sd <- c( old[1], old[2], old[3] )
		R  <- matrix( c( 1, old[4], old[5], old[4], 1, old[6], old[5], old[6], 1 ), 3, 3 )
		accept_sd  <- rep( 0, 3 )
		accept_rho <- rep( 0, 3 )

		#- helper: full log-posterior kernel for given (sd, R)
		logpost <- function( sd, R ) {
			SIGMA <- diag(sd) %*% R %*% diag(sd)
			-(sum_np/2) * base::determinant( SIGMA, logarithm=TRUE )$modulus[1] -
				0.5 * sum( RES * base::solve( SIGMA ) ) +
				sum( bayes_dhalft( sd, scale=scale, nu=nu, log=TRUE ) ) +
				bayes_dlkjcorr( R, eta=eta, log=TRUE )
		}

		#- (A) update each SD individually via log-RW (Jacobian sign corrected):
		for ( i in 1:3 ) {

			sd_c <- sd
			logsd_p   <- rnorm( 1, mean=log( sd[i] ), sd=tau_sd )
			sd_p      <- sd
			sd_p[i]   <- exp( logsd_p )

			llratio <- ( logpost( sd_p, R ) - logpost( sd_c, R ) ) +
				( logsd_p - log( sd_c[i] ) )    # correct Jacobian: +log(new) - log(old)

			if ( log( runif(1) ) < llratio ) {
				sd <- sd_p
				accept_sd[i] <- 1
			}
		}

		#- (B) update each correlation individually within its exact PD-feasible range:
		pairs  <- rbind( c(1,2), c(1,3), c(2,3) )   # (j,k) to update
		others <- rbind( c(1,3,2,3), c(1,2,2,3), c(1,2,1,3) )  
		# others[m,]: indices (jm1,jm2,km1,km2) s.t. a=R[jm1,jm2], b=R[km1,km2] -- see below

		for ( m in 1:3 ) {

			jk <- pairs[m,]; j <- jk[1]; k <- jk[2]
			idx3 <- setdiff( 1:3, jk )           # the "third" index

			a <- R[j, idx3]   # rho_{j, third}
			b <- R[k, idx3]   # rho_{k, third}

			#- exact PD-feasible interval for rho_{jk} given a, b fixed:
			root <- sqrt( (1 - a^2) * (1 - b^2) )
			lo <- a*b - root
			hi <- a*b + root

			rho_c <- R[j,k]
			u_c   <- ( rho_c - lo ) / ( hi - lo )
			x_c   <- log( u_c / (1 - u_c) )             # logit

			x_p   <- rnorm( 1, mean=x_c, sd=tau_rho )
			u_p   <- 1 / ( 1 + exp(-x_p) )               # sigmoid
			rho_p <- lo + (hi - lo) * u_p

			R_p <- R
			R_p[j,k] <- R_p[k,j] <- rho_p

			jac <- ( log( u_p * (1 - u_p) ) ) - ( log( u_c * (1 - u_c) ) )

			llratio <- ( logpost( sd, R_p ) - logpost( sd, R ) ) + jac

			if ( log( runif(1) ) < llratio ) {
				R <- R_p
				accept_rho[m] <- 1
			}
		}

		#- finalize:
		SIGMA_P <- diag(sd) %*% R %*% diag(sd)
		update  <- c( sd, R[ rbind( c(1,2), c(1,3), c(2,3) ) ] )
		accept  <- c( accept_sd, accept_rho )

	}

	#- output:
	return( list( SIGMA_P=SIGMA_P, update_p=update, accept_p=accept ) )

}

tsrm_bayes_update_sigmad <- function( dyads=NULL, args_list=NULL, sum_nd=NULL ) 
{

	#- make transformation matrix Q:
	c2 <- 1 / sqrt(2)
	Q <- matrix(0, 6, 6)
	Q[1, c(1, 2)] <- c(c2,  c2)   # s_ab
	Q[2, c(3, 4)] <- c(c2,  c2)   # s_ac
	Q[3, c(5, 6)] <- c(c2,  c2)   # s_bc
	Q[4, c(1, 2)] <- c(c2, -c2)   # d_ab
	Q[5, c(3, 4)] <- c(c2, -c2)   # d_ac
	Q[6, c(5, 6)] <- c(c2, -c2)   # d_bc

	#- transform random effects:
	rs <- matrix( do.call( "c", dyads ), ncol=6, byrow=TRUE )
	us <- rs%*%t(Q)

	#- compute residual sum of square matrix ( sum u*u^T )
	RES_s <- t( us[,c(1:3)] )%*%( us[,c(1:3)] )
	RES_d <- t( us[,c(4:6)] )%*%( us[,c(4:6)] )

	#- get hyperparms:
	nud <- args_list$mcmc_args$IW.nud
	Sd1 <- args_list$mcmc_args$IW.Sd[[1]] 
	Sd2 <- args_list$mcmc_args$IW.Sd[[2]]  

	#- generate draw:
	PHI_D <- matrix(0, 6, 6)
	PHI_D[1:3,1:3] <- MCMCpack::riwish( nud + sum_nd, Sd1 + RES_s )
	PHI_D[4:6,4:6] <- MCMCpack::riwish( nud + sum_nd, Sd2 + RES_d )

	#- compute SIGMA_D:
	SIGMA_D <- t(Q)%*%PHI_D%*%Q
	R_D     <- cov2cor( SIGMA_D )
	update  <- c( sqrt( diag( SIGMA_D )[c(1,3,5)] ), 
		R_D[rbind( c(1,2),c(3,4),c(5,6),c(1,3),c(1,6),c(3,5),c(1,4),c(1,5),c(3,6) )] )
	accept  <- 1

	#- output:
	return( list( SIGMA_D=SIGMA_D, update_d=update, accept=accept ) )

}

tsrm_bayes_update_sigmat <- function(  data_list=NULL, args_list=NULL, 
	BETA=NULL, person=NULL, dyads=NULL ) 
{

	#- make transformation matrix Q:
	s6 <- sqrt(6); s3 <- sqrt(3)
	Q  <- matrix(c(
		  1/s6,  1/s6,  1/s6,      1/s6,      1/s6,      1/s6,      
		  1/s6, -1/s6, -1/s6,      1/s6,      1/s6,     -1/s6,      
		  1/s3, -1/s3,  1/(2*s3), -1/(2*s3), -1/(2*s3),  1/(2*s3),  
		  0,     0,     0.5,      -0.5,       0.5,      -0.5,       
		  0,     0,     0.5,       0.5,      -0.5,      -0.5,       
		  1/s3,  1/s3, -1/(2*s3), -1/(2*s3), -1/(2*s3), -1/(2*s3)   
	), nrow = 6, byrow = TRUE)

	#- we need to compute the residuals e_ijk
	groupinfo <- data_list[["groupinfo"]]
  	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	Zt <- data_list[["Zt"]]
	nt <- data_list[["nt"]]

	#- get no. groups:
	ngroups <- nrow( groupinfo )

	#- iterate:
	residuals <- vector("list", ngroups )
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

		#- compute residual:
		resi <- as.vector( yi - ( Xi%*%BETA + Zpi%*%peffi + Zdi%*%deffi ) )

		residuals[[ng]] <- as.vector( t( Zti) %*% resi ) # for the order
			
	}

	#- transform residuals to a matrix:
	es <- matrix( do.call( "c", residuals ), ncol=6, byrow=TRUE )
	us <- es%*%t(Q)
	sum_nt <- nrow( us )

	#- compute residual sum of square matrix ( sum u*u^T )
	RES_1 <- sum( us[,1]^2 )
	RES_2 <- sum( us[,2]^2 )
	#RES_A <- t( us[,c(3:4)] )%*%( us[,c(3:4)] ) + t( us[,c(5:6)] )%*%( us[,c(5:6)] )
	RES_A <- t( us[,c(3,5)] )%*%( us[,c(3,5)] ) + t( us[,c(4,6)] )%*%( us[,c(4,6)] )

	#- get hyperparms:
	at  <- args_list$mcmc_args$IG.at
	bt  <- args_list$mcmc_args$IG.bt
	nut <- args_list$mcmc_args$IW.nut 
	St  <- args_list$mcmc_args$IW.St 
		
	#- generate draw:
	PHI_T <- matrix(0, 6, 6)
	PHI_T[1,1] <- 1/rgamma(1, nut + 0.5*sum_nt, bt + 0.5*RES_1 )
	PHI_T[2,2] <- 1/rgamma(1, nut + 0.5*sum_nt, bt + 0.5*RES_2 )
	# PHI_T[3:4,3:4] <- PHI_T[5:6,5:6] <- MCMCpack::riwish( nut + 2*sum_nt, St + RES_A )
	A_draw <- MCMCpack::riwish( nut + 2*sum_nt, St + RES_A )
	PHI_T[3,3] <- A_draw[1,1]; PHI_T[5,5] <- A_draw[2,2]
	PHI_T[3,5] <- PHI_T[5,3] <- A_draw[1,2]
	PHI_T[4,4] <- A_draw[1,1]; PHI_T[6,6] <- A_draw[2,2]
	PHI_T[4,6] <- PHI_T[6,4] <- A_draw[1,2]

	#- compute SIGMA_T:
	SIGMA_T <- t(Q)%*%PHI_T%*%Q
	R_T     <- cov2cor( SIGMA_T )
	update  <- c( sqrt( SIGMA_T[1,1] ), R_T[rbind( c(1,2),c(1,6),c(1,3),c(1,4) )] )
	accept  <- 1

	#- output:
	return( list( SIGMA_T=SIGMA_T, update_t=update, accept=accept ) )

}
