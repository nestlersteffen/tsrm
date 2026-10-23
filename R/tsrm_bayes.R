
#---- a wrapper function for the MCMC sampler for the TSRM parms

tsrm_bayes <- function( parm_table=NULL, parm_list=NULL, data_list=NULL, 
	names_list=NULL, tsrm_data=NULL, args_list=NULL, nchain=NULL, biter=NULL, 
	burnin=NULL, verbose=TRUE, cl=NULL )
{
	
	#- insert parm
	parm_list <- tsrm_include_free_parameters( parm=parm_table$starts, parm_list=parm_list, 
		parm_table=parm_table )

	#- get covariance matrices:
	parm_list <- get_sigmas( parm_list=parm_list )

	#- get initial estimates of person and dyad effects:
	initial_effects <- tsrm_bayes_initialize( tsrm_data=tsrm_data, data_list=data_list, 
		names_list=names_list, args_list=args_list )
	person <- initial_effects[["person"]]
	dyads  <- initial_effects[["dyads"]]

	#- update hyperparameters:
	args_list <- tsrm_bayes_hyperparms( args_list=args_list, data_list=data_list, 
		person=person, dyads=dyads )

	#- sampling:
    if ( nchain == 1 ) {

    	result <- tsrm_bayes_chain( inits=parm_table$starts, parm_list=parm_list, 
    		data_list=data_list, args_list=args_list, person=person, dyads=dyads, 
    		biter=biter, verbose=verbose )
    	result <- list( result )

    } else {

    	if (is.null(cl)) {
	        cl <- parallel::makeCluster(nchain)
	        doParallel::registerDoParallel(cl)
	        on.exit(parallel::stopCluster(cl), add = TRUE)

	        parallel::clusterEvalQ(cl, {
	            RhpcBLASctl::blas_set_num_threads(1)
	            RhpcBLASctl::omp_set_num_threads(1)
	        })

	    } else {
	        doParallel::registerDoParallel(cl)
	    }

        #- run chains in parallel:
	    nc <- NULL
	    result <- foreach::foreach( nc = seq(1, nchain, 1 ), .packages = "tsrm" ) %dopar% {
        
        	tsrm_bayes_chain( inits=parm_table$starts, parm_list=parm_list, 
	    		data_list=data_list, args_list=args_list, person=person, dyads=dyads, 
	    		biter=biter, verbose=FALSE )
    	} 
        
    } 

    #- make a cube with the sampled thetas and with deleted burnin:
    parms 		<- array( 0, dim = c( biter, nrow( parm_table ), nchain ) )
    parms_final <- array( 0, dim = c( biter-burnin, nrow( parm_table ), nchain ) )
    for ( nc in 1:nchain ) {
    	parms[,,nc]       <- result[[nc]]$theta
        parms_final[,,nc] <- parms[-c(1:burnin),,nc]
    }

    #- compute convergence diagnostics:
    iact <- matrix( 0, nrow = nchain, ncol = nrow( parm_table ) )
    ess  <- matrix( 0, nrow = nchain, ncol = nrow( parm_table ) )
    for ( nc in 1:nchain ) {
    	iact[nc,]  <- apply( parms_final[,,nc], 2, bayes_iact )
    	ess[nc,]   <- nrow( parms_final[,,nc] )*(1/iact[nc,] )
    }
    iact  <- colMeans( iact )
    ess.m <- colMeans( ess )
    ess.o <- colSums( ess )

    rhat <- NULL
    if ( nchain > 1 ) {
        rhat <- bayes_rhat( parms_final )
    } 

    #- compute acceptance probs
    p_accept <- p_accept_post <- vector("list", nchain )
    for ( nc in 1:nchain ) {
    	p_accept[[nc]] <- colSums( result[[nc]]$accept )/ biter
    	p_accept_post[[nc]] <- colSums( result[[nc]]$accept[-(1:burnin),,drop=FALSE] )/ (biter-burnin)
    }
    p_accept <- Reduce("+", p_accept)/nchain
    p_accept_post <- Reduce("+", p_accept_post)/nchain

	#- get final estimates:
	#parms        	<- do.call( rbind, lapply( 1:nchain, function(nc) parms[,,nc] ) )
	parms_final  	<- do.call( rbind, lapply( 1:nchain, function(nc) parms_final[,,nc] ) )
	parm_table$p.M  <- apply( parms_final, 2, mean )
	parm_table$p.Md <- apply( parms_final, 2, quantile, prob = 0.50 )
	parm_table$p.SD <- apply( parms_final, 2, sd )
	parm_table$lCI  <- apply( parms_final, 2, quantile, prob = 0.025 )
	parm_table$uCI  <- apply( parms_final, 2, quantile, prob = 0.975 )
	parm_table$iact <- iact
	parm_table$ess.m <- ess.m
	parm_table$ess.o <- ess.o
	parm_table$rhat  <- if ( nchain > 1 ) rhat
	
	out <- list( parm_table=parm_table, p_accept=p_accept, p_accept_post=p_accept_post, 
		chains=parms, args_list=args_list )
	return( out )
}