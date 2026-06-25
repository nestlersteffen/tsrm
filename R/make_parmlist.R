
#---- this function generates the parameter matrices per group:

make_parmlist <- function( parm_table=NULL, names_list=NULL, args_list=NULL, model=c("srm","tsrm") )
{
	#- get no. of variables:
    no_vars <- names_list$no_var
	
	#- make BETA based on number of BETA parms in parm_table
	idx 		 <- which( parm_table$type == "BETA" )
	lgt_mu_preds <- length( idx )
	BETA         <- matrix( 0.0, nrow = lgt_mu_preds, ncol = 1 )

	#- make the different covraiance, correlations etc. matrices
	if ( model == "srm" ) {
		SD_G  	<- diag( 0, no_vars*1 )
		SD_P  	<- diag( 1, no_vars*2 )
		SD_D  	<- diag( 1, no_vars*2 )
		SD_T  	<- matrix( 0, nrow=1, ncol=0 )
		SIGMA_G <- RHO_G <- diag( 1, no_vars*1 )
		SIGMA_P <- RHO_P <- diag( 1, no_vars*2 )
		SIGMA_D <- RHO_D <- diag( 1, no_vars*2 )
		SIGMA_T <- RHO_T <- matrix( 0, nrow=1, ncol=0 )
	} else if ( model == "tsrm" ) {
		SD_G  	<- matrix( 0, nrow=1, ncol=0 )
		SD_P  	<- diag( 0, no_vars*3 )
		SD_D  	<- diag( 0, no_vars*6 )
		SD_T  	<- diag( 0, no_vars*6 )
		SIGMA_G <- RHO_G <- matrix( 0, nrow=1, ncol=0 )
		SIGMA_P <- RHO_P <- diag( 1, no_vars*3 )
		SIGMA_D <- RHO_D <- diag( 1, no_vars*6 )
		SIGMA_T <- RHO_T <- diag( 1, no_vars*6 )
	} else if ( model == "htsrm" ) {
		SD_G  	<- matrix( 0, nrow=1, ncol=0 )
		SD_P  	<- diag( 0, 5 )
		SD_D  	<- diag( 0, 8 )
		SD_T  	<- diag( 0, 6 )
		SIGMA_G <- RHO_G <- matrix( 0, nrow=1, ncol=0 )
		SIGMA_P <- RHO_P <- diag( 1, 3 )
		SIGMA_D <- RHO_D <- diag( 1, 8 )
		SIGMA_T <- RHO_T <- diag( 1, 6 )
	}
	
	parm_list <- list( BETA=BETA, SD_G=SD_G, SD_P=SD_P, SD_D=SD_D, SD_T=SD_T, 
		RHO_G=RHO_G, RHO_P=RHO_P, RHO_D=RHO_D, RHO_T=RHO_T, 
		SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, SIGMA_T=SIGMA_T )

	#- consider the random group:
	return( parm_list )
} 