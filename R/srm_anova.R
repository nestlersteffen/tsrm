
srm_anova_singlegroup <- function( group_data=NULL, np=NULL, names_list=NULL, with_info=FALSE ) 
{

	#- get relevant variables
	p_var   <- names_list$p_var
	d_var   <- names_list$d_var
	d_type  <- names_list$d_var_type
	outcome <- names_list$outcome
	no_var  <- names_list$no_var 

	#- ---------------------------------------------
	#-  Step 1: compute different (centered) means

	person <- matrix(0,nrow=np,ncol=2*no_var)
	dyads  <- matrix(0,nrow=np*(np-1),ncol=2*no_var)
	for ( nv in seq( no_var ) ) {
		
		#- make a measure-specific dataset
		tmp_data <- group_data[ group_data[,"measure"] == nv, ]
		
		#- compute grand mean
		grand <- mean( tmp_data[,outcome] )
		
		#- compute person-level effects
		person[,2*nv - 1 ] <- anova_means( tmp_data[,outcome], tmp_data[,p_var[1]] ) - grand
		person[,2*nv] <- anova_means( tmp_data[,outcome], tmp_data[,p_var[2]] ) - grand

		#- dyad - level effects:
		idx <- which( tmp_data[,d_type[1]] == 1 )
		ab_ij <- anova_means( tmp_data[idx,outcome], tmp_data[idx,d_var[1]] ) - grand
		ab_ji <- anova_means( tmp_data[-idx,outcome], tmp_data[-idx,d_var[1]] ) - grand
		dyads[,2*nv-1] <- c( ab_ij, ab_ji )
		dyads[,2*nv]   <- c( ab_ji, ab_ij )

	}

	#- ---------------------------------
	#-  Step 2: compute sum of squares

	#- constant factors
	f1 <- np - 1
	
	#- person-level variance and covariance terms:
	psumsq <- f1*crossprod(person)
	
	#- dyad-level variance and covariance terms:
	dsumsq <- crossprod(dyads)

	#- -------------------------------------------
	#-  Step 3: get S and obtain final estimates
	S <- srm_anova_buildS( N=np )

	if ( no_var == 1 ) {
		sumsq <- numeric(6) 
		sumsq <- c( psumsq[2,2], psumsq[1,1], dsumsq[1,1], psumsq[1,2], psumsq[1,2], dsumsq[1,2])
		gparm <- S%*%sumsq
		gparm <- gparm[c(2,1,4,3,6)]
	} else {
		sumsq <- numeric(18) 
		sumsq <- c( psumsq[2,2], psumsq[1,1], dsumsq[1,1], psumsq[1,2], psumsq[1,2], dsumsq[1,2],
					psumsq[3,3], psumsq[4,4], dsumsq[3,3], psumsq[3,4], psumsq[3,4], dsumsq[3,4],
					psumsq[2,3], psumsq[1,4], dsumsq[1,4], psumsq[1,3], psumsq[2,4], dsumsq[1,3] )
		gparm <- as.vector(S %*% matrix(sumsq, nrow = 6, ncol = 3))
		gparm <- gparm[c(2,1,7,8,4,16,14,13,17,3,9,6,18,15,12)]
	}

	#- return ...
	if ( with_info ) {
		return( list( gparm=gparm, person=person, dyads=dyads ) )
	} else {
		return( gparm )
	}
	
}