
tsrm_anova_singlegroup <- function( group_data=NULL, np=NULL, names_list=NULL ) 
{

	#- get relevant variables
	p_var   <- names_list$p_var
	d_var   <- names_list$d_var
	d_type  <- names_list$d_var_type
	t_var   <- names_list$t_var
	t_type  <- names_list$t_var_type
	outcome <- names_list$outcome

	#- make empty vector:
	sumsq <- numeric(33)

	#- ---------------------------------
	#-  Step 1: compute different (centered) means
	grand   <- mean( group_data[,outcome] )
	
	#- person - level:
	actor   <- anova_means( group_data[,outcome], group_data[,p_var[1]] ) - grand
	partner <- anova_means( group_data[,outcome], group_data[,p_var[2]] ) - grand
	judge   <- anova_means( group_data[,outcome], group_data[,p_var[3]] ) - grand
	
	#- dyad - level:	
	idx <- which( group_data[,d_type[1]] == 1 )
	ab_ij <- anova_means( group_data[idx,outcome], 
		group_data[idx,d_var[1]] ) - grand
	ab_ji <- anova_means( group_data[-idx,outcome], 
		group_data[-idx,d_var[1]] ) - grand
	
	idx <- which( group_data[,d_type[2]] == 1 )
	ac_ij <- anova_means( group_data[idx,outcome], 
		group_data[idx,d_var[2]] ) - grand
	ac_ji <- anova_means( group_data[-idx,outcome], 
		group_data[-idx,d_var[2]] ) - grand
	
	idx <- which( group_data[,d_type[3]] == 1 )
	bc_ij <- anova_means( group_data[idx,outcome], 
		group_data[idx,d_var[3]] ) - grand
	bc_ji <- anova_means( group_data[-idx,outcome], 
		group_data[-idx,d_var[3]] ) - grand
	
	#- triad - level 
	abc_kij <- anova_means( group_data[group_data[,t_type] == 1,outcome], 
		group_data[group_data[,t_type] == 1,t_var] ) - grand
	abc_kji <- anova_means( group_data[group_data[,t_type] == 2,outcome], 
		group_data[group_data[,t_type] == 2,t_var] ) - grand
	abc_jik <- anova_means( group_data[group_data[,t_type] == 3,outcome], 
		group_data[group_data[,t_type] == 3,t_var] ) - grand
	abc_jki <- anova_means( group_data[group_data[,t_type] == 4,outcome], 
		group_data[group_data[,t_type] == 4,t_var] ) - grand
	abc_ijk <- anova_means( group_data[group_data[,t_type] == 5,outcome], 
		group_data[group_data[,t_type] == 5,t_var] ) - grand
	abc_ikj <- anova_means( group_data[group_data[,t_type] == 6,outcome], 
		group_data[group_data[,t_type] == 6,t_var] ) - grand

	#- ---------------------------------
	#-  Step 2: compute sum of squares

	#- constant factors
	f1 <- np - 1
	f2 <- np - 2

	# person-level (variances and covariances):
	sumsq[1] <- f1*f2*sum( judge^2 )   # judge variance
  	sumsq[2] <- f1*f2*sum( actor^2 )   # perceiver variance
  	sumsq[3] <- f1*f2*sum( partner^2 ) # target variance
  	sumsq[4] <- sumsq[6] <- f1*f2*sum( judge*actor )   # judge-perceiver covariance
  	sumsq[5] <- sumsq[8] <- f1*f2*sum( judge*partner ) # judge-target covariance
  	sumsq[7] <- sumsq[9] <- f1*f2*sum( actor*partner ) # perceiver-target covariance
  
  	# dyad-level (variances and covariances):
  	sumsq[10] <- f2*sum( c( ac_ij, ac_ji )^2 ) #judge x perceiver variance
	sumsq[11] <- f2*sum( c( bc_ij, bc_ji )^2 ) #judge x target variance
	sumsq[12] <- f2*sum( c( ab_ij, ab_ji )^2 ) #perceiver x target variance
	sumsq[13] <- f2*sum( c( ac_ij, ac_ji )*c( ac_ji, ac_ij ) ) # judge x perceiver - perceiver x judge covariance
  	sumsq[14] <- f2*sum( c( bc_ij, bc_ji )*c( bc_ji, bc_ij ) ) # judge x target - target x judge covariance
	sumsq[15] <- f2*sum( c( ab_ij, ab_ji )*c( ab_ji, ab_ij ) ) # perceiver x target - target x perceiver covariance
  	sumsq[16] <- sumsq[18] <- f2*sum( c( ac_ij, ac_ji )*c( bc_ij, bc_ji ) ) # judge x perceiver - judge x target covariance
   	sumsq[17] <- sumsq[20] <- f2*sum( c( ac_ij, ac_ji )*c( ab_ij, ab_ji ) ) # judge x perceiver - perceiver x target covariance
  	sumsq[19] <- sumsq[21] <- f2*sum( c( bc_ij, bc_ji )*c( ab_ji, ab_ij ) ) # judge x target - target x perceiver covariance
	sumsq[22] <- sumsq[24] <- f2*sum( c( ac_ij, ac_ji )*c( bc_ji, bc_ij ) ) # judge x perceiver - target x judge covariance
  	sumsq[23] <- sumsq[26] <- f2*sum( c( ab_ij, ab_ji )*c( ac_ji, ac_ij ) ) # perceiver x target - perceiver x judge covariance
  	sumsq[25] <- sumsq[27] <- f2*sum( c( bc_ij, bc_ji )*c( ab_ij, ab_ji ) ) # judge x target - target x perceiver covariance

  	# triad_level (variances and covariances ):
  	sumsq[28] <- sum( ( group_data[,outcome] - grand )^2 )
  	sumsq[29] <- sum( # within judge
  		c( abc_ijk, abc_jik, abc_ikj, abc_kij, abc_jki, abc_kji )* 
		c( abc_jik, abc_ijk, abc_kij, abc_ikj, abc_kji, abc_jki ) )
	sumsq[31] <- sum( # within targets
		c( abc_ijk, abc_jik, abc_ikj, abc_kij, abc_jki, abc_kji )* 
		c( abc_ikj, abc_jki, abc_ijk, abc_kji, abc_jik, abc_kij ) )
	sumsq[30] <- sum( # within perceivers
		c( abc_ijk, abc_jik, abc_ikj, abc_kij, abc_jki, abc_kji )* 
		c( abc_kji, abc_kij, abc_jki, abc_jik, abc_ikj, abc_ijk ) )
	if ( np != 5 ) {
		sumsq[32] <- sumsq[33] <- sum( # pure within
			c( abc_ijk, abc_jik, abc_ikj, abc_kij, abc_jki, abc_kji )* 
			c( (abc_kij+abc_jki)/2,(abc_ikj+abc_kji)/2,(abc_jik+abc_kji)/2, 
			   (abc_ijk+abc_jki)/2,(abc_ijk+abc_kij)/2,(abc_jik+abc_ikj)/2) ) 
	}
	
	#- -------------------------------------------
	#-  Step 3: get S and obtain final estimates
	S    <- tsrm_anova_buildS( N=np )
	invS <- base::solve(S)

	return( invS%*%sumsq )

}

tsrm_anova_pool <- function( parms=NULL, group_of_5=NULL, parm_table=NULL ) 
{
	#- pool the results:
	ngroups   <- nrow( parms)
	parm_mean <- colMeans( parms )
   	parm_ses  <- apply( parms, 2, sd )/sqrt( ngroups )
   	if ( any( group_of_5 ) ) {
   		no_group_of_5 <- sum( group_of_5 )
   		parm_mean[32:33] <- colMeans( parms[-group_of_5,32:33] )
   		parm_ses[32:33]  <- apply( parms[-group_of_5,32:33], 2, sd )/sqrt( ngroups - no_group_of_5 )
   	}

   	#- get final estimates:
   	parm_order <- c(1,2,3,10,11,12,28,4,5,7,13,14,15,16,19,23,22,17,25,31,30,29,32)
   	est <- parm_mean[parm_order]
   	ses <- parm_ses[parm_order]

   	#- make a table
   	tab <- data.frame(Est=est, Std.Error=ses)
   	return( list( tab=tab, parm_mean=parm_mean, parm_ses=parm_ses ) )

}