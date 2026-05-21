
#' Fit a Triadic Social Relations Model to Generalized Round-Robin Data
#'
#' @param formula A formula for the fixed effects
#' @param g_var A vector with an identifier for round-robin groups
#' @param p_var A vector with identifiers for actor, partner, and judge
#' @param data A data frame
#' @param control A list of control arguments
#' @param debug A control args for debugging, for internal use only
#' @export

tsrm <- function( formula = NULL, p_var = NULL, g_var = NULL, data = NULL, 
	control=make_args_list(), debug = FALSE )
{

	mm <- match.call()
	
	#--- step 0: make args_list:
	missCtrl <- missing( control )
	if ( !missCtrl && !inherits( control, "make_args_list" ) ) {
  		if(!is.list(control)) { stop("'control' has to be a list.") }
  		args_list <- do.call( make_args_list, control )
 	} else {
 		args_list <- control
 	}
 	
	#--- step 1: get formulas and variables names
	mm_formula <- mm[["formula"]]
	if (is.null(mm_formula)) {
  		stop("formula parameter is missing.")
	}

	if (!any(grepl("$", mm_formula, fixed=TRUE))) {
  		formula <- eval(mm_formula, data, enclos=sys.frame(sys.parent()))
	} else {
  		formula <- eval(mm_formula, envir = parent.frame())
	}

	if (inherits(formula, "formula")) {
  		# univariate case - make a list
  		formulas <- list(formula)
	} else if (is.list(formula)) {
  		# bi-and multivariate case - formula is a list
  		formulas <- formula
	} else {
  		stop("formula must be a formula or a list of formulas.")
	}

	names_y <- c()
	names_X <- list()

	for (i in seq_along(formulas)) {
  		if (!inherits(formulas[[i]], "formula")) {
    		stop(paste("i", "th formula is not valid."))
  		}
  		names_y[i]   <- all.vars(formulas[[i]])[1]
  		names_X[[i]] <- attr(stats::terms(formulas[[i]]), "term.labels")
	}
	
	#--- step 2: make tsrm_data_frame 
	tmp <- tsrm_make_dataframe( names_y=names_y, names_X=names_X, p_var=p_var,
   	 	g_var=g_var, formulas=formulas, data=data )
  	tsrm_data   <- tmp$data
  	names_list  <- tmp$names_list

  	#--- step 3: make data_list
	data_list <- make_datalist( data=tsrm_data, names_list=names_list, 
	 	args_list=args_list, model="tsrm" ) 
  	
  	#--- step 4: make parm_table:
	parm_table <- make_parmtable( data_list=data_list, names_list=names_list, 
	 	args_list=args_list, model="tsrm" )

	#--- step 5: make parm_list:
	parm_list  <- make_parmlist( parm_table=parm_table, names_list=names_list, 
		args_list=args_list, model="tsrm" )

	#--- step 6: add start values ( optional )
	parm_table <- add_starts( parm_table=parm_table, data_list=data_list, 
	 	args_list=args_list )

	if ( debug ) {
		result <- list( parm_table=parm_table, parm_list=parm_list, data_list=data_list, 
			args_list=args_list, names_list=names_list, data=tsrm_data )
		return( result )
	}

	#--- step 7: fit the model
	result <- fit_model( parm_table=parm_table, parm_list=parm_list, data_list=data_list, 
		args_list=args_list, model="tsrm" ) 

	#- some warnings: 
  	if ( !result$converged ) {
  	   warning("The algorithm did not converge.")
  	}
  	if ( result$warning_vcov ) {
  	   warning("Hessian is not invertible. The covariance matrix of the 
  	   	estimates could not be computed.")
  	}
  
	#- finally...
	result <- c( result, model="tsrm", list( names_list=names_list, names_y=names_y ) )
	class( result ) <- "tsrm"
  	return( result )

}