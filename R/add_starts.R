
#---- this function generates some starting values

add_starts <- function( data=NULL, parm_table=NULL, parm_list=NULL, data_list=NULL, 
	args_list=NULL, names_list=NULL, model=c("srm","tsrm","htsrm") ) 
{

	#- define start value-, lower- and upper-bound-vector 
  	starts <- numeric( dim( parm_table )[1] )
  	low    <- rep( -Inf, dim( parm_table )[1] )
  	up     <- rep(  Inf, dim( parm_table )[1] )

  	#- should we use anova-estimates?
  	if ( args_list$starts == "anova" ) {
 		anova_table <- anova( data=data, names_list=names_list, parm_table=parm_table, 
 			model=model)
 		anova_table <- anova_transform_estimates( parm_table=anova_table, parm_list=parm_list, 
 			names_list=names_list, model=model, standardization=FALSE)
	}

	#- start values for BETA:
	X 	    <- data_list$X
	y 	    <- data_list$y
	fit_lm  <- stats::lm(y~0+X)
	betas   <- coef( fit_lm )
	if ( args_list$fixed_group ) {
		betas <- c( rep( betas[1], nrow( data_list$groupinfo) ), betas[-1] )
	} 
	idx <- which( parm_table$type == "BETA" )
	starts[idx] <- betas

	#- start values for the variance terms:
	#  SIGMA_G/SD_G:
	idx <- which( parm_table$type %in% c("SIGMA_G","SD_G")  & 
	   			  parm_table$pos1 == parm_table$pos2 )
	starts[idx] <- 0.1#log( 1 )

	#- start values for variance in SIGMA_P/SD_P
	idx <- which( parm_table$type %in% c("SIGMA_P","SD_P") & 
				  parm_table$pos1 == parm_table$pos2 )
	starts[idx] <- if ( args_list$starts == "anova" ) anova_table[idx,"anova_transformed"] else 1

	#- start values for variance in SIGMA_D/SD_D
	idx <- which( parm_table$type %in% c("SIGMA_D","SD_D") & 
				  parm_table$pos1 == parm_table$pos2 )
	starts[idx] <- if ( args_list$starts == "anova" ) anova_table[idx,"anova_transformed"] else 1

	#- start values for variance in SIGMA_T/SD_T
	idx <- which( parm_table$type %in% c("SIGMA_T","SD_T")  & 
	   			  parm_table$pos1 == parm_table$pos2 )
	starts[idx] <- if ( args_list$starts == "anova" ) anova_table[idx,"anova_transformed"] else 1

	#- output:
	parm_table$starts <- starts
	parm_table$low    <- low
  	parm_table$up     <- up
	return( parm_table )

}