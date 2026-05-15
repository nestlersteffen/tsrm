
#---- this function generates the design matrices per group:

srm_make_parmtable <- function( data_list = NULL, names_list = NULL, 
	args_list = NULL ) 
{
	
	#- get no. of variables:
	no_var <- names_list$no_var

	#- parm_table for the univariate SRM with predictors:

	if ( no_var == 1L ) {

		#- everything for BETA:
		mu_preds 	 <- names_list[["mu_preds"]]
		if ( args_list$fixed_group ) {
			# get number of groups:
			ng <- nrow( data_list$groupinfo )
			# number of predictors ( minus 1 because mu_preds contains the intercept )
			lgt_mu_preds <- ng + length( do.call( "c", mu_preds ) ) - 1
		} else {
			lgt_mu_preds <- length( do.call( "c", mu_preds ) ) 
		}
		BETA_type  <- rep("BETA", lgt_mu_preds )
		BETA_pos1  <- seq(1,lgt_mu_preds,1)
		BETA_pos2  <- rep(1,lgt_mu_preds)
		BETA_ntype <- rep(0,lgt_mu_preds)

		#- everything for SIGMA_P/S_P and RHO_P:
		S_P_type  <- rep("SD_P", 2 )
		S_P_pos1  <- c(1,2)
		S_P_pos2  <- c(1,2)
		S_P_ntype <- rep(1,2)
		R_P_type  <- rep("RHO_P", 1 )
		R_P_pos1  <- c(1)
		R_P_pos2  <- c(2)
		R_P_ntype <- rep(2,1)

		#- everything for SIGMA_D/S_D and RHO_D:
		S_D_type  <- rep("SD_D", 1 )
		S_D_pos1  <- c(1)
		S_D_pos2  <- c(1)
		S_D_ntype <- rep(3,1)
		R_D_type  <- rep("RHO_D", 1 )
		R_D_pos1  <- c(1)
		R_D_pos2  <- c(2)
		R_D_ntype <- rep(4,1)

	} else {

		#- everything for BETA:
		mu_preds 	 <- names_list[["mu_preds"]]
		lgt_mu_preds <- length( do.call( "c", mu_preds ) )
		BETA_type    <- rep("BETA", lgt_mu_preds )
		BETA_pos1    <- seq(1,lgt_mu_preds,1)
		BETA_pos2    <- rep(1,lgt_mu_preds)
		BETA_ntype   <- rep(0,lgt_mu_preds)

		#- everything for SIGMA_P/S_P and RHO_P:
		S_P_type <- rep("SD_P", 4 )
		S_P_pos1 <- c(1,2,3,4)
		S_P_pos2 <- c(1,2,3,4)
		S_P_ntype <- rep(1,4)
		R_P_type <- rep("RHO_P", 6 )
		R_P_pos1 <- c(1,1,1,2,2,3)
		R_P_pos2 <- c(2,3,4,3,4,4)
		R_P_ntype <- rep(2,6)

		#- everything for SIGMA_D/S_D and RHO_D:
		S_D_type  <- rep("SD_D", 2 )
		S_D_pos1  <- c(1,3)
		S_D_pos2  <- c(1,3)
		S_D_ntype <- rep(3,2)
		R_D_type  <- rep("RHO_D", 4 )
		R_D_pos1  <- c(1,1,1,3)
		R_D_pos2  <- c(2,3,4,4)
		R_D_ntype <- rep(4,4)

	}

	parm_table <- data.frame( 
	   type  = c( BETA_type, S_P_type, R_P_type, S_D_type, R_D_type ),
	   pos1  = c( BETA_pos1, S_P_pos1, R_P_pos1, S_D_pos1, R_D_pos1 ),
	   pos2  = c( BETA_pos2, S_P_pos2, R_P_pos2, S_D_pos2, R_D_pos2 ),
	   ntype = c( BETA_ntype, S_P_ntype, R_P_ntype, S_D_ntype, R_D_ntype ) 
	)

	if ( args_list$random_group ) {
		
		if ( no_var == 1L ) {

			parm_table_group <- data.frame( 
		 		type  = c( "SD_G" ),
	     		pos1  = c( 1 ),
	     		pos2  = c( 1 ),
	     		ntype = c( 5 )
	    	)
		
		} else {
			parm_table_group <- data.frame( 
		 		type  = c( rep( "SD_G", 2 ), "RHO_G" ),
	     		pos1  = c( 1, 2, 1 ),
	     		pos2  = c( 1, 2, 2 ),
	     		ntype = c( rep( 5, 2 ), 6 )
	    	)
		}

		parm_table <- rbind( parm_table, parm_table_group )

	}

	# some final things:
	parm_table$index <- seq(1, nrow(parm_table), 1 )
	parm_table$type  <- as.character( parm_table$type )
	return( parm_table )
}