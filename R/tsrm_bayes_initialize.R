
#------- function to initialize some things for TSRM sampler

tsrm_bayes_initialize <- function( tsrm_data = NULL, data_list = NULL, names_list = NULL, 
	args_list = NULL )
{

	#- get variable names:
	g_var    <- names_list[["g_var"]]
	p_var   <- names_list$p_var
	d_var   <- names_list$d_var
	d_type  <- names_list$d_var_type
	t_var   <- names_list$t_var
	t_type  <- names_list$t_var_type
	outcome <- names_list$outcome
	no_var  <- names_list$no_var 
	
	#- get infos for the groups:
	groups    <- unique( tsrm_data[,g_var] )
	ngroups   <- base::length( groups)
	nv        <- data_list[["nv"]]

	#- we generate person effect estimates per group
	person <- vector("list", ngroups )
	dyads  <- vector("list", ngroups )

	#- iterate across groups to add columns if necessary:
	for ( ng in seq( ngroups ) ) {

		#- get group-specific data-frame:
		tmp_data <- tsrm_data[tsrm_data[,g_var] == groups[ng] & tsrm_data[,"measure",] == nv,]

		#- compute grand mean
		grand <- mean( tmp_data[,outcome] )

		#- person - level:
		actor 	<- anova_means( tmp_data[,outcome], tmp_data[,p_var[1]] ) - grand
		partner <- anova_means( tmp_data[,outcome], tmp_data[,p_var[2]] ) - grand
		judge 	<- anova_means( tmp_data[,outcome], tmp_data[,p_var[3]] ) - grand
	
		#- dyad - level effects:
		idx <- which( tmp_data[,d_type[1]] == 1 )
		ab_ij <- anova_means( tmp_data[idx,outcome], tmp_data[idx,d_var[1]] ) - grand
		ab_ji <- anova_means( tmp_data[-idx,outcome], tmp_data[-idx,d_var[1]] ) - grand
	
		idx <- which( tmp_data[,d_type[2]] == 1 )
		ac_ij <- anova_means( tmp_data[idx,outcome], tmp_data[idx,d_var[2]] ) - grand
		ac_ji <- anova_means( tmp_data[-idx,outcome], tmp_data[-idx,d_var[2]] ) - grand
	
		idx <- which( tmp_data[,d_type[3]] == 1 )
		bc_ij <- anova_means( tmp_data[idx,outcome], tmp_data[idx,d_var[3]] ) - grand
		bc_ji <- anova_means( tmp_data[-idx,outcome], tmp_data[-idx,d_var[3]] ) - grand

		#- save the elements in a list:
		person[[ng]] <- c( t( cbind( actor, partner, judge ) ) )
		dyads[[ng]] <- c( t( cbind( ab_ij, ab_ji, ac_ij, ac_ji, bc_ij, bc_ji ) ) )

	}

	#-
	return( list( person = person, dyads = dyads ) )

}