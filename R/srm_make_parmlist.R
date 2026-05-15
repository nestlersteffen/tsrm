
#---- this function generates the design matrices per group:

srm_make_parmlist <- function( parm_table = NULL, names_list = NULL, 
	args_list = NULL ) 
{
	#- get some names:
	no_vars <- names_list$no_var
	
	#- get number of BETA parms in parm_table
	idx 		 <- which( parm_table$type == "BETA" )
	lgt_mu_preds <- length( idx )
	#if ( no_vars != 1L ) {
	#	lgt_mu_preds <- length( do.call( "c", names_list[["mu_preds"]] ) )
	#} else {
	#	lgt_mu_preds <- length( names_list[["mu_preds"]] )
	#}
	
	#- a default parm_list:
	BETA  <- matrix( 0.0, nrow = lgt_mu_preds, ncol = 1 )
	SD_G  <- diag( 0.0, no_vars*1 )
	SD_P  <- diag( 1.0, no_vars*2 )
	SD_D  <- diag( 1.0, no_vars*2 )
	RHO_G <- diag( 1.0, no_vars*1 )
	RHO_P <- diag( 1.0, no_vars*2 )
	RHO_D <- diag( 1.0, no_vars*2 )
	parm_list <- list( BETA = BETA, SD_G = SD_G, SD_P = SD_P, SD_D = SD_D, 
		RHO_G = RHO_G, RHO_P = RHO_P, RHO_D = RHO_D, SIGMA_G = NULL, SIGMA_P = NULL, 
		SIGMA_D = NULL )

	#- result
	return( parm_list )
}

