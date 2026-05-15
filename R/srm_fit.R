
#---- this function contains routines to fit the SRM:

srm_fit <- function( parm_table = NULL, parm_list = NULL, data_list = NULL,
	args_list = NULL )
{

	#- get the start values
	starts    <- parm_table$starts
	
	#- the algorithm has not converged...
	converged    <- FALSE
	warning_vcov <- TRUE

	#- ---------------------------
	#-     fit the model

	if ( args_list$method == "c++" ) {
		
		#- make the optimization function:
		obj <- tsrm_MakeOptFun( srm_compute_gradient, parm_list=parm_list, 
			parm_table=parm_table, data_list=data_list, 
			args_list=args_list, both=TRUE ) 

		#- get variance estimates:
		fit <- suppressWarnings( tryCatch( 
			stats::nlminb( start=starts, objective=obj$fn, gradient=obj$gr, 
				control = args_list$nlminb_ctrl ),
		    error   = function(e) { NULL }
		) )
		
	} else if ( args_list$method == "tmb" ) {

		#- make tmb-data_list:
		tdata <- list( y = data_list$y, X = data_list$X, Zp = data_list$Zp, 
			Zd = data_list$Zd, np = data_list$np, nd = data_list$nd, 
			groupinfo = data_list$groupinfo, method = args_list$with_reml, 
			no_var = data_list$nv )
		
		#- make tmb-parm_list
		idx_BETA <- which( parm_table$type == "BETA")
		idx_sdP  <- which( parm_table$type == "SD_P")
		idx_rhoP <- which( parm_table$type == "RHO_P")
		idx_sdD  <- which( parm_table$type == "SD_D")
		idx_rhoD <- which( parm_table$type == "RHO_D")
		tparm <- list( BETA=starts[idx_BETA], sdP=starts[idx_sdP], 
			rhoP=starts[idx_rhoP], sdD=starts[idx_sdD], rhoD=starts[idx_rhoD] )

		#- some adds in case of random groups:
		if ( args_list$random_group ) {
			
			#- add to tdata_list:
		 	tdata <- c( tdata, Zg = data_list$Zg, model = "srm_groups_tmb" )
		 	
		 	#- add to tparm
		 	idx_sdG  <- which( parm_table$type == "SD_G")
		 	if ( data_list$nv == 1 ) {
		 		rhoG <- c(0)
		 	} else if ( data_list$nv == 2 ) {
		 		idx_rhoG <- which( parm_table$type == "RHO_G")
		 		rhoG     <- starts[idx_rhoG]
		 	}
		 	tparm <- c( tparm, sdG = starts[idx_sdG], rhoG = rhoG )

	 	} else { 
		 	
		 	tdata <- c( tdata, model = "srm_tmb" )
		
		}

		#- make tmb object:
		obj <- TMB::MakeADFun( data = tdata, parameters = tparm, 
			DLL="tsrm_TMBExports", silent=TRUE )

	 	#- fit model:
		fit <- suppressWarnings( tryCatch( 
			stats::nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, 
				control=args_list$nlminb_ctrl ),
		    error   = function(e) { NULL }
		) )
	}

	#- error handling:
	if (  is.null( fit ) || inherits( fit, "try-error" ) ) {
	    parm_new <- starts
	    break
	} else {
		converged <- TRUE
		parm_new  <- fit$par
	    ll        <- -1*fit$objective
	    parm_list <- srm_include_free_parameters( parm=parm_new, 
	    	parm_list=parm_list, parm_table=parm_table )
	}  

	#- --------------------------------------------------

	if ( converged ) {
		
		#- get standard errors:
		parm_table <- srm_get_ses( parm=parm_new, parm_table=parm_table, 
			parm_list=parm_list, data_list=data_list, args_list=args_list ) 

		#- deviance and aic
		dev <- -2*ll
		aic <- dev + 2*length( parm_new )
		warning_vcov   <- FALSE
		warning_nlminb <- fit$convergence

	} else { 
		ll <- dev <- aic <- vcov <- NULL
	}

	#- save and return: 
  	res <- list( parm = parm_new, parm_table = parm_table, parm_list = parm_list,
    	data_list = data_list, args_list = args_list, deviance = dev, ll = ll, 
    	aic = aic, vcov = vcov, converged = converged, warning_vcov = warning_vcov,
    	warning_nlminb=warning_nlminb )            
  	return( res )
}