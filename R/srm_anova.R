
srm_anova_singlegroup <- function( group_data=NULL, np=NULL, names_list=NULL ) 
{

	#- get relevant variables
	p_var   <- names_list$p_var
	d_var   <- names_list$d_var
	d_type  <- names_list$d_var_type
	outcome <- names_list$outcome

	#- make empty vector:
	sumsq <- numeric(6)

	#- ---------------------------------
	#-  Step 1: compute different (centered) means
	grand   <- mean( group_data[,outcome] )
	
	#- person - level:
	actor   <- anova_means( group_data[,outcome], group_data[,p_var[1]] ) - grand
	partner <- anova_means( group_data[,outcome], group_data[,p_var[2]] ) - grand
	
	#- dyad - level:	
	idx <- which( group_data[,d_type[1]] == 1 )
	ab_ij <- anova_means( group_data[idx,outcome], 
		group_data[idx,d_var[1]] ) - grand
	ab_ji <- anova_means( group_data[-idx,outcome], 
		group_data[-idx,d_var[1]] ) - grand
	
	#- ---------------------------------
	#-  Step 2: compute sum of squares

	#- constant factors
	f1 <- np - 1
	
	# person-level (variances and covariances):
	sumsq[2] <- f1*sum( actor^2 )   # perceiver variance
  	sumsq[1] <- f1*sum( partner^2 ) # target variance
  	sumsq[4] <- sumsq[5] <- f1*sum( actor*partner ) # perceiver-target covariance
  
  	# dyad-level (variances and covariances):
  	sumsq[3] <- sum( c( ab_ij, ab_ji )^2 ) # perceiver x target variance
	sumsq[6] <- sum( c( ab_ij, ab_ji )*c( ab_ji, ab_ij ) ) # perceiver x target - target x perceiver covariance
	
	#- -------------------------------------------
	#-  Step 3: get S and obtain final estimates
	S <- srm_anova_buildS( N=np )
	return( S%*%sumsq )

}

srm_anova_pool <- function( parms=NULL, parm_table=NULL, group_of_5=NULL ) 
{
	#- pool the results:
	ngroups   <- nrow( parms)
	parm_mean <- colMeans( parms )
   	
	#- get the standard errors:
	# if ( args_list$anova_ses = "LB" && no_var == 1 ) {
	# 	next
	# } else {
		parm_ses  <- apply( parms, 2, sd )/sqrt( ngroups )
	# }
   	
   	#- get final estimates:
   	parm_order <- c(2,1,4,3,6)
   	est <- parm_mean[parm_order]
   	ses <- parm_ses[parm_order]

   	#- make a table
   	tab <- data.frame(Est=est, Std.Error=ses)
   	return( list( tab=tab, parm_mean=parm_mean, parm_ses=parm_ses ) )

}