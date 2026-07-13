

anova <- function( data=NULL, names_list=NULL, parm_table=NULL, model=c("srm","tsrm","htsrm"))
{

	#- get model-specific functions.
	model <- match.arg( model )
	if ( model == "srm" ) {
	 	anova_singlegroup <- srm_anova_singlegroup
	 	anova_pool 		  <- srm_anova_pool
	 	if ( names_list$no_var == 1 ) no_parms <- 6 else no_parms <- 18
	} else if ( model == "tsrm" ) {
		anova_singlegroup <- tsrm_anova_singlegroup
		anova_pool 		  <- tsrm_anova_pool
		if ( names_list$no_var == 1 ) no_parms <- 33 else no_parms <- 79
	} else if ( model == "htsrm" ) {
		anova_singlegroup <- htsrm_anova_singlegroup
		anova_pool 		  <- htsrm_anova_pool
		no_parms          <- 43
	} else {
		stop("False model class defined (anova).")
	}

	#- get groupinfo matrix:

	#- get names of relevant variables:
	g_var   <- names_list$g_var
	p_var   <- names_list$p_var
	
	#- get no. groups:
	group_ids  <- unique( data[,g_var] )
	ngroups    <- length( group_ids )
	parms      <- matrix( 0, nrow=ngroups, ncol=no_parms)
	group_of_5 <- rep(FALSE, ngroups)

	#- compute group-specific parameter estimates:
	for ( ng in seq( ngroups ) ) {
			
		#- get group-specific dataframe:
		group_data <- data[ data[ ,g_var] == group_ids[ ng ], ]

		#- check whether group contains five persons only
		np <- length( unique( unlist( group_data[, p_var] ) ) )
		if ( np == 5 ) group_of_5[ng] <- TRUE

		#- comoute parms:
		parms[ng,] <- anova_singlegroup( group_data=group_data, np=np, names_list=names_list)

   	}

   	#- pool the group-specific estimates:
   	result <- anova_pool( parms=parms, group_of_5=group_of_5, parm_table=parm_table ) 
 	return( result )
}