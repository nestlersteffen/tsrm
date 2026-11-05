
#---- this function contains routines to fit the SRM:

fit_model <- function( parm_table=NULL, parm_list=NULL, data_list=NULL, args_list=NULL, 
	model=c("srm","tsrm","htsrm") )
{

	#- get model and the correct function:
	model <- match.arg( model )
	include_free_parms <- switch( model,
        "srm"  = srm_include_free_parameters,
        "tsrm" = tsrm_include_free_parameters,
        "htsrm" = htsrm_include_free_parameters,
        stop( paste( "The srm, the tsrm, or the htsrm can be estimated (fit)." ) )
    )

	#- get the start values
	starts    <- parm_table$starts
	
	#- the algorithm has not converged...
	converged    <- FALSE
	warning_vcov <- TRUE

	#- ---------------------------
	#-     fit the model

	if ( !args_list$constrained ) {

		#print("ohne constraint")
		
		#- make the optimization function:
		obj <- make_optfct( compute_gradfct, parm_list=parm_list, 
			parm_table=parm_table, data_list=data_list, 
			args_list=args_list, model=model, both=TRUE ) 

		#- get variance estimates:
		fit <- suppressWarnings( tryCatch( 
			stats::nlminb( start=starts, objective=obj$fn, gradient=obj$gr, 
				control = args_list$nlminb_ctrl ),
		    error   = function(e) { NULL }
		) )
		
	} else {

		#print("mit constraint")

		#- make tmb-data_list:
		tdata <- list( y = data_list$y, X = data_list$X, Zp = data_list$Zp, 
			Zd = data_list$Zd, np = data_list$np, nd = data_list$nd, 
			groupinfo = data_list$groupinfo, method = args_list$with_reml )
		
		#- make tmb-parm_list
		idx_BETA <- which( parm_table$type == "BETA")
		idx_sdP  <- which( parm_table$type == "SD_P")
		idx_rhoP <- which( parm_table$type == "RHO_P")
		idx_sdD  <- which( parm_table$type == "SD_D")
		idx_rhoD <- which( parm_table$type == "RHO_D")
		tparm <- list( BETA=starts[idx_BETA], vecsdP=log(starts[idx_sdP]), vecrhoP=starts[idx_rhoP] )

		#- some adds in case of random groups:
		if ( model == "srm" ) {

			tparm <- append( tparm, list( vecsdD=log(starts[idx_sdD]), vecrhoD=starts[idx_rhoD] ) ) 

			if ( args_list$random_group ) {
			
				#- add to tdata_list:
		 		tdata <- append( tdata, list( Zg = data_list$Zg ), 2 )
		 		tdata <- c( tdata, model = "srm_groups_tmb" )
		 	
		 		#- add to tparm
		 		idx_sdG  <- which( parm_table$type == "SD_G")
		 		if ( data_list$nv == 1 ) {
		 			rhoG <- c(0)
		 		} else if ( data_list$nv == 2 ) {
		 			idx_rhoG <- which( parm_table$type == "RHO_G")
		 			rhoG     <- starts[idx_rhoG]
		 		}
		 		tparm <- c( tparm, vecsdG = log( starts[idx_sdG] ) )

	 		} else { 
		 	
		 		tdata <- c( tdata, model = "srm_tmb" )
		
			}
		
		} else {

			#- add to data_list:
			tdata <- append( tdata, list( Zt = data_list$Zt ), 4 )
			tdata <- append( tdata, list( nt = data_list$nt ), 7 )
			tdata <- c( tdata, list( version = 1, model = "tsrm_tmb" ) )

			#- add to tparm:
			idx_sdT  <- which( parm_table$type == "SD_T")
			idx_rhoT <- which( parm_table$type == "RHO_T")
			tparm    <- append( tparm, list( 
				vecsdD=rep( log( starts[idx_sdD] ), times=2 ), vecrhoD=starts[idx_rhoD],
				vecsdT=starts[idx_sdT], vecrhoT=starts[idx_rhoT] ) )

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
	    warning("Something bad happened during optimization.")
	} else {
		converged <- TRUE
		parm_new  <- fit$par
	    ll        <- -1*fit$objective
	}  

	#- --------------------------------------------------

	if ( converged ) {
		
		#- get standard errors:
		if ( !args_list$constrained ) {
			#print("ohne constraint")
			parm_table <- get_ses( parm=parm_new, parm_table=parm_table, parm_list=parm_list, 
				data_list=data_list, args_list=args_list, model=model ) 
		} else {
			#print("mit constraint")
			tmb_results   <- tryCatch({  TMB::sdreport( obj ) },
				error = function(e) NULL )
			if ( !is.null( tmb_results ) ) {
				tmb_estimates <- unname( tmb_results$value )
				tmb_stderrors <- tmb_results$sd
			} else {
				tmb_estimates <- tmb_stderrors <- rep( NA, length( parm_new ) )
				warning("Something bad happened in TMB's computations.")
			}
			parm_table$est <- tmb_estimates
			parm_table$se  <- tmb_stderrors
		}

		#- update parm_list:
		parm_list <- include_free_parms( parm=parm_table$est, parm_list=parm_list, 
			parm_table=parm_table )

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