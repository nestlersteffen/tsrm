
tsrm_anova_singlegroup <- function( group_data=NULL, np=NULL, names_list=NULL, with_info=FALSE ) 
{

	#- get relevant variables
	p_var   <- names_list$p_var
	d_var   <- names_list$d_var
	d_type  <- names_list$d_var_type
	t_var   <- names_list$t_var
	t_type  <- names_list$t_var_type
	outcome <- names_list$outcome
	no_var  <- names_list$no_var 

	#- ---------------------------------
	#-  Step 1: compute different (centered) means
	
	person <- matrix(0,nrow=np,ncol=3*no_var)
	dyads  <- matrix(0,nrow=np*(np-1),ncol=6*no_var)
	triads <- matrix(0,nrow=np*(np-1)*(np-2),ncol=5*no_var)

	for ( nv in seq( no_var ) ) {
		
		#- make a measure-specific dataset
		tmp_data <- group_data[ group_data[,"measure"] == nv, ]
		
		#- compute grand mean
		grand <- mean( tmp_data[,outcome] )

		#- person - level:
		person[,1 + 3*(nv-1) ] <- anova_means( tmp_data[,outcome], tmp_data[,p_var[1]] ) - grand
		person[,2 + 3*(nv-1) ] <- anova_means( tmp_data[,outcome], tmp_data[,p_var[2]] ) - grand
		person[,3 + 3*(nv-1) ] <- anova_means( tmp_data[,outcome], tmp_data[,p_var[3]] ) - grand
	
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

		dyads[,1 + 6*(nv-1)] <- c( ab_ij, ab_ji )
		dyads[,2 + 6*(nv-1)] <- c( ab_ji, ab_ij )
		dyads[,3 + 6*(nv-1)] <- c( ac_ij, ac_ji )
		dyads[,4 + 6*(nv-1)] <- c( ac_ji, ac_ij )
		dyads[,5 + 6*(nv-1)] <- c( bc_ij, bc_ji )
		dyads[,6 + 6*(nv-1)] <- c( bc_ji, bc_ij )

		#- triad - level 
		abc_kij <- anova_means( tmp_data[tmp_data[,t_type] == 1,outcome], 
			tmp_data[tmp_data[,t_type] == 1,t_var] ) - grand
		abc_kji <- anova_means( tmp_data[tmp_data[,t_type] == 2,outcome], 
			tmp_data[tmp_data[,t_type] == 2,t_var] ) - grand
		abc_jik <- anova_means( tmp_data[tmp_data[,t_type] == 3,outcome], 
			tmp_data[tmp_data[,t_type] == 3,t_var] ) - grand
		abc_jki <- anova_means( tmp_data[tmp_data[,t_type] == 4,outcome], 
			tmp_data[tmp_data[,t_type] == 4,t_var] ) - grand
		abc_ijk <- anova_means( tmp_data[tmp_data[,t_type] == 5,outcome], 
			tmp_data[tmp_data[,t_type] == 5,t_var] ) - grand
		abc_ikj <- anova_means( tmp_data[tmp_data[,t_type] == 6,outcome], 
			tmp_data[tmp_data[,t_type] == 6,t_var] ) - grand

		triads[,1 + 5*(nv-1)] <- c( abc_ijk, abc_jik, abc_ikj, abc_kij, abc_jki, abc_kji )
		triads[,2 + 5*(nv-1)] <- c( abc_jik, abc_ijk, abc_kij, abc_ikj, abc_kji, abc_jki )
		triads[,3 + 5*(nv-1)] <- c( abc_kji, abc_kij, abc_jki, abc_jik, abc_ikj, abc_ijk )
		triads[,4 + 5*(nv-1)] <- c( abc_ikj, abc_jki, abc_ijk, abc_kji, abc_jik, abc_kij )
		triads[,5 + 5*(nv-1)] <- c( (abc_kij+abc_jki)/2,(abc_ikj+abc_kji)/2,(abc_jik+abc_kji)/2, 
									(abc_ijk+abc_jki)/2,(abc_ijk+abc_kij)/2,(abc_jik+abc_ikj)/2 )
	}

	#- ---------------------------------
	#-  Step 2: compute sum of squares

	#- constant factors
	f1 <- np - 1
	f2 <- np - 2

	#- person-level variance and covariance terms:
	psumsq <- f1*f2*crossprod(person)

	#- dyad-level variance and covariance terms:
	dsumsq <- f2*crossprod(dyads)

	#- dyad-level variance and covariance terms:
	tsumsq <- crossprod(triads)

	#- -------------------------------------------
	#-  Step 3: get S and obtain final estimates
	
	S     <- tsrm_anova_buildS( N=np )
	invS  <- base::solve(S)

	if ( no_var == 1 ) {
		sumsq <- c( psumsq[3,3], psumsq[1,1], psumsq[2,2], psumsq[1,3], psumsq[2,3], 
	 				psumsq[1,3], psumsq[1,2], psumsq[2,3], psumsq[1,2],
					dsumsq[3,3], dsumsq[5,5], dsumsq[1,1], dsumsq[3,4], dsumsq[5,6], dsumsq[1,2],
					dsumsq[3,5], dsumsq[3,1], dsumsq[3,5], dsumsq[5,2], dsumsq[3,1], dsumsq[5,2],
					dsumsq[3,6], dsumsq[1,4], dsumsq[3,6], dsumsq[1,5], dsumsq[1,4], dsumsq[1,5],
					
					tsumsq[1,1], tsumsq[1,2], tsumsq[1,3], tsumsq[1,4], tsumsq[1,5], tsumsq[1,5] )
		gparm <- invS%*%sumsq
		gparm <- gparm[c(2,3,1,7,6,5,12,10,11,15,13,14,17,21,16,23,25,24,28,31,29,30,32)]
	} else {
		#- split the matrices
		psumsq1  <- psumsq[1:3,1:3]
		psumsq2  <- psumsq[4:6,4:6]
		psumsq12 <- psumsq[1:3,4:6]
		dsumsq1  <- dsumsq[1:6,1:6]
		dsumsq2  <- dsumsq[7:12,7:12]
		dsumsq12 <- dsumsq[1:6,7:12]
		tsumsq1  <- tsumsq[1:5,1:5]
		tsumsq2  <- tsumsq[6:10,6:10]
		tsumsq12 <- tsumsq[1:5,6:10]
		#- make vectors
		sumsq <- c( # Variable 1
					psumsq1[3,3], psumsq1[1,1], psumsq1[2,2], psumsq1[1,3], psumsq1[2,3], 
	 				psumsq1[1,3], psumsq1[1,2], psumsq1[2,3], psumsq1[1,2],
					dsumsq1[3,3], dsumsq1[5,5], dsumsq1[1,1], dsumsq1[3,4], dsumsq1[5,6], dsumsq1[1,2],
					dsumsq1[3,5], dsumsq1[3,1], dsumsq1[3,5], dsumsq1[5,2], dsumsq1[3,1], dsumsq1[5,2],
					dsumsq1[3,6], dsumsq1[1,4], dsumsq1[3,6], dsumsq1[1,5], dsumsq1[1,4], dsumsq1[1,5],
					tsumsq1[1,1], tsumsq1[1,2], tsumsq1[1,3], tsumsq1[1,4], tsumsq1[1,5], tsumsq1[1,5],
					# Variable 2
					psumsq2[3,3], psumsq2[1,1], psumsq2[2,2], psumsq2[1,3], psumsq2[2,3], 
	 				psumsq2[1,3], psumsq2[1,2], psumsq2[2,3], psumsq2[1,2],
					dsumsq2[3,3], dsumsq2[5,5], dsumsq2[1,1], dsumsq2[3,4], dsumsq2[5,6], dsumsq2[1,2],
					dsumsq2[3,5], dsumsq2[3,1], dsumsq2[3,5], dsumsq2[5,2], dsumsq2[3,1], dsumsq2[5,2],
					dsumsq2[3,6], dsumsq2[1,4], dsumsq2[3,6], dsumsq2[1,5], dsumsq2[1,4], dsumsq2[1,5],
					tsumsq2[1,1], tsumsq2[1,2], tsumsq2[1,3], tsumsq2[1,4], tsumsq2[1,5], tsumsq2[1,5],
					# Variable 3 (ab 67)
					psumsq12[3,3], psumsq12[1,1], psumsq12[2,2], psumsq12[3,1], psumsq12[3,2], 
	 				psumsq12[1,3], psumsq12[1,2], psumsq12[2,3], psumsq12[2,1],
	 				dsumsq12[3,3], # 76 --
	 				dsumsq12[5,5], # 77 --
	 				dsumsq12[1,1], # 78 --
	 				dsumsq12[3,4], # 79 --
	 				dsumsq12[5,6], # 80 --
	 				dsumsq12[1,2], # 81 --
					dsumsq12[3,5], # 82 --
					dsumsq12[3,1], # 83 --
					dsumsq12[5,3], # 84 --
					dsumsq12[5,2], # 85 -- 
					###
					dsumsq12[1,4], # 86 -- # dsumsq12[3,1], 
					dsumsq12[2,5], # 87 
					dsumsq12[3,6], # 88 --
					dsumsq12[4,1], # 89 (5-10)
					dsumsq12[6,3], # 90 (3-8)
					dsumsq12[5,1], # 91 --
					dsumsq12[1,3], # 92 -- # dsumsq12[1,4], 
					dsumsq12[1,5], # 93 --
					tsumsq12[1,1], tsumsq12[1,2], tsumsq12[1,3], tsumsq12[1,4], tsumsq12[2,4], tsumsq12[2,3] )
		gparm <- as.vector(invS %*% matrix(sumsq, nrow = 33, ncol = 3))
		gparm <- gparm[c(2,3,1,35,36,34,
			7,6,5,40,39,38, 68,73,72,75,69,74,70,71,67,
			12,10,11,45,43,44,
			15,13,14,17,21,16,23,25,24, 48,46,47,50,54,49,56,58,57,
			78,81,92,86,93,87, 
			83,90,76,79,82,88,
			91,85,84,89,77,80,
			28,61,31,29,30,32,64,62,63,65,
			94,96,98,95,97,99 )]

	}
	
	#- return ...
	if ( with_info ) {
		return( list( gparm=gparm, person=person, dyads=dyads ) )
	} else {
		return( gparm )
	}

}