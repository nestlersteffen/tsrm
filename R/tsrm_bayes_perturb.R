
#---- function to generate chain-specific starting values for the TSRM Bayes sampler

tsrm_bayes_perturb <- function( parm_list=NULL, person=NULL, dyads=NULL, args_list=NULL, 
	effects=c( "prior", "keep" ) )
{

	#- a function used to enforce a symmetric matrix:
	sym <- function( M ) ( M + t(M) ) / 2

	#- checks and defaults:
	effects  <- match.arg( effects )
	sd_log   <- args_list$mcmc_args$sd_log
	eta      <- args_list$mcmc_args$eta
	nu       <- args_list$mcmc_args$nu
	
	#- SIGMA_P: perturb the SDs, draw a new correlation matrix from LKJ
	sd_old  <- sqrt( diag( parm_list[["SIGMA_P"]] ) )
	sd_new  <- sd_old * exp( rnorm( 3, 0, sd_log ) )
	R_new   <- matrix( as.numeric( bayes_rlkjcorr( 1, 3, eta=eta ) ), 3, 3 )
	new_SIGMA_P <- sym( diag( sd_new, 3 ) %*% R_new %*% diag( sd_new, 3 ) )

	#- SIGMA_D: perturb the two 3 x 3 blocks of Phi_D (same Q as in tsrm_bayes_update_sigmad)
	c2 <- 1 / sqrt(2)
	Q_D <- matrix( 0, 6, 6 )
	Q_D[1, c(1, 2)] <- c( c2,  c2 )   # s_ab
	Q_D[2, c(3, 4)] <- c( c2,  c2 )   # s_ac
	Q_D[3, c(5, 6)] <- c( c2,  c2 )   # s_bc
	Q_D[4, c(1, 2)] <- c( c2, -c2 )   # d_ab
	Q_D[5, c(3, 4)] <- c( c2, -c2 )   # d_ac
	Q_D[6, c(5, 6)] <- c( c2, -c2 )   # d_bc

	PHI_D_old <- Q_D %*% parm_list[["SIGMA_D"]] %*% t( Q_D )
	PHI_D     <- matrix( 0, 6, 6 )
	PHI_D[1:3,1:3] <- MCMCpack::riwish( nu, (nu-4) * sym( PHI_D_old[1:3,1:3] ) )
	PHI_D[4:6,4:6] <- MCMCpack::riwish( nu, (nu-4) * sym( PHI_D_old[4:6,4:6] ) )
	new_SIGMA_D    <- sym( t( Q_D ) %*% PHI_D %*% Q_D )

	#- SIGMA_T: perturb Phi_T (same Q as in tsrm_bayes_update_sigmat)

	s6 <- sqrt(6); s3 <- sqrt(3)
	Q_T <- matrix( c(
		  1/s6,  1/s6,  1/s6,      1/s6,      1/s6,      1/s6,
		  1/s6, -1/s6, -1/s6,      1/s6,      1/s6,     -1/s6,
		  1/s3, -1/s3,  1/(2*s3), -1/(2*s3), -1/(2*s3),  1/(2*s3),
		  0,     0,     0.5,      -0.5,       0.5,      -0.5,
		  0,     0,     0.5,       0.5,      -0.5,      -0.5,
		  1/s3,  1/s3, -1/(2*s3), -1/(2*s3), -1/(2*s3), -1/(2*s3)
	), nrow = 6, byrow = TRUE )

	PHI_T_old <- Q_T %*% parm_list[["SIGMA_T"]] %*% t( Q_T )

	# the two scalar blocks: inverse gamma with expectation phi_old (shape nu/2, rate (nu-2)*phi/2)
	phi1 <- 1 / rgamma( 1, shape = nu/2, rate = (nu-2) * PHI_T_old[1,1] / 2 )
	phi2 <- 1 / rgamma( 1, shape = nu/2, rate = (nu-2) * PHI_T_old[2,2] / 2 )

	# the 2 x 2 block that appears twice (positions (3,5) and (4,6)); average both copies:
	A_old <- ( PHI_T_old[c(3,5),c(3,5)] + PHI_T_old[c(4,6),c(4,6)] ) / 2
	A_new <- MCMCpack::riwish( nu, (nu-3) * sym( A_old ) )

	PHI_T <- matrix( 0, 6, 6 )
	PHI_T[1,1] <- phi1
	PHI_T[2,2] <- phi2
	PHI_T[3,3] <- A_new[1,1]; PHI_T[5,5] <- A_new[2,2]
	PHI_T[3,5] <- PHI_T[5,3] <- A_new[1,2]
	PHI_T[4,4] <- A_new[1,1]; PHI_T[6,6] <- A_new[2,2]
	PHI_T[4,6] <- PHI_T[6,4] <- A_new[1,2]
	new_SIGMA_T <- sym( t( Q_T ) %*% PHI_T %*% Q_T )

	#- perturb effects  (or keep them)

	new_person <- person
	new_dyads  <- dyads

	if ( effects == "prior" ) {

		for ( ng in seq_along( person ) ) {
			n_p <- length( person[[ng]] ) / 3
			M   <- mvtnorm::rmvnorm( n_p, sigma = new_SIGMA_P )   # n_p x 3
			new_person[[ng]] <- c( t( M ) )                        # (A_1,B_1,C_1,A_2,...)
		}

		for ( ng in seq_along( dyads ) ) {
			n_d <- length( dyads[[ng]] ) / 6
			M   <- mvtnorm::rmvnorm( n_d, sigma = new_SIGMA_D )   # n_d x 6
			new_dyads[[ng]] <- c( t( M ) )                         # 6 entries per dyad
		}

	}

	#- output:
	parm_list[["SIGMA_P"]] <- new_SIGMA_P
	parm_list[["SIGMA_D"]] <- new_SIGMA_D
	parm_list[["SIGMA_T"]] <- new_SIGMA_T
	return( list( parm_list = parm_list, person = new_person, dyads = new_dyads ) )

}