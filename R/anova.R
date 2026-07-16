

anova <- function( data=NULL, names_list=NULL, parm_table=NULL, with_ses=TRUE,
	model=c("srm","tsrm","htsrm"))
{

	#- get model-specific functions.
	model <- match.arg( model )
	if ( model == "srm" ) {
	 	anova_singlegroup <- srm_anova_singlegroup
	 	if ( names_list$no_var == 1 ) no_parms <- 5 else no_parms <- 16
	} else if ( model == "tsrm" ) {
		anova_singlegroup <- tsrm_anova_singlegroup
		if ( names_list$no_var == 1 ) no_parms <- 23 else no_parms <- 79
	} else if ( model == "htsrm" ) {
		anova_singlegroup <- htsrm_anova_singlegroup
		no_parms <- 40
	} else {
		stop("False model class defined (anova).")
	}

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

   	#- get the average of the group-specific estimates:
   	parm_mean <- colMeans( parms )
   	if ( model != "srm" & any( group_of_5 ) ) {
   		parm_mean[no_parms] <- colMeans( parms[-group_of_5,no_parms] )
   	}

   	#- compute the standard errors of the average estimates:
   	parm_ses <- rep( NA, no_parms )
   	if ( with_ses & ngroups > 1 ) {
   		parm_ses <- apply( parms, 2, sd )/sqrt( ngroups )
   		if ( model != "srm" & any( group_of_5 ) ) {
   			parm_ses[no_parms] <- apply( parms[-group_of_5,no_parms], 2, sd )/sqrt( ngroups - no_group_of_5 )
   		}
   	}

   	#- add to parm_table:
   	idx <- which( parm_table$type=="BETA" )
   	parm_table$anova_est <- c( rep(NA,length(idx)), parm_mean )
   	parm_table$anova_se  <- c( rep(NA,length(idx)), parm_ses )

   	#- make output object:
   	# result <- data.frame(Est=parm_mean, Std.Error=parm_ses) 
 	return( parm_table )
}