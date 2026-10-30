
#---- function that implements a chain of a MCMC sampler for the TSRM parms

tsrm_bayes_chain <- function( parm_table=NULL, parm_list=NULL, data_list=NULL, args_list=NULL, 
	person=NULL, dyads=NULL, biter=NULL, verbose=NULL )
{
	
	#- collect parameter matrices:
	BETA    <- parm_list[["BETA"]]
	SIGMA_P <- parm_list[["SIGMA_P"]]
	SIGMA_D <- parm_list[["SIGMA_D"]]
	SIGMA_T <- parm_list[["SIGMA_T"]]
	sum_np  <- sum( data_list$groupinfo[,1] )
	sum_nd  <- sum( data_list$groupinfo[,2] )
	
	#- matrices to save the draws and the acceptance rates:
	theta  <- matrix( 0, nrow=biter, length( parm_table$starts ) )
	accept <- vector( "list", biter)

	#- the first value for update_p
	update_p  <- c( sqrt( diag( SIGMA_P ) ), rep(0,3) ) 
	
	#- do some sampling...
	for ( xx in 1:biter ) {

		if ( verbose & xx%%50 == 0 ) {
        	print( paste0( "Iteration: ", xx ) )
    	}

		#- update BETA:
		# Rprof("prof.out", interval = 0.005, line.profiling = TRUE)
		tmp_BETA <- tsrm_bayes_update_beta( data_list=data_list, args_list=args_list, 
			person=person, dyads=dyads, SIGMA_T=SIGMA_T )
		# Rprof(NULL)
		BETA 		 <- tmp_BETA$BETA
		OBJ_invVepsi <- tmp_BETA$OBJ_invVepsi # a list
		
		#- update p and r:
		effects <- tsrm_bayes_update_effects( data_list=data_list, args_list=args_list, 
			person=person, dyads=dyads, BETA=BETA, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, 
			OBJ_invVepsi=OBJ_invVepsi )
		person <- effects[["person"]]
		dyads  <- effects[["dyads"]]
		
		#- update SIGMA_P:
		tmp_sigmap <- tsrm_bayes_update_sigmap( person=person, args_list=args_list, sum_np=sum_np, 
			old=update_p ) 
		SIGMA_P    <- tmp_sigmap$SIGMA_P
		update_p   <- tmp_sigmap$update_p
		accept_p   <- tmp_sigmap$accept_p
		
		#- update SIGMA_D:
		tmp_sigmad <- tsrm_bayes_update_sigmad( dyads=dyads, args_list=args_list, sum_nd=sum_nd) 
		SIGMA_D    <- tmp_sigmad$SIGMA_D
		update_d   <- tmp_sigmad$update_d
		
	    #- update SIGMA_T:
	    tmp_sigmat <- tsrm_bayes_update_sigmat( data_list=data_list, args_list=args_list, 
			BETA=BETA, person=person, dyads=dyads ) 
	    SIGMA_T    <- tmp_sigmat$SIGMA_T
		update_t   <- tmp_sigmat$update_t
		
		#- save all things:
		theta[xx,]   <- c( BETA, update_p, update_d, update_t )
		accept[[xx]] <- accept_p

	}

	#- return result
	result <- list( theta=theta, accept=do.call("rbind", accept ) )
	return( result )
}