
#---- this function generates the design matrices per group:

srm_make_dataframe <- function( names_y = NULL, names_X = NULL, p_var = NULL,
  	d_var = NULL, g_var = NULL, formulas = NULL, data = NULL ) 
{

	data_names <- colnames( data )

	#- some checks...
	if ( !all( p_var %in% data_names ) | is.null( p_var ) ) {
		stop("Something is wrong with the person-level identifiers.")
	}
	if ( !is.null( d_var ) ) {
		if ( !( d_var %in% data_names ) ) {
		stop("Something is wrong with the dyad-level identifier.")
		}
	}
	if ( !( all( names_y %in% data_names ) ) ) {
		stop("Something is wrong with the outcome variable.")
	}

	#- we build a dataframe:
	srm_data <- data[,c(p_var,d_var,names_y)]
	
	#- add a group-variable if necessary:
	if ( is.null( g_var ) ) {
	  	srm_data$group <- 1
	  	g_var <- "group"
	} else {
		srm_data <- cbind( srm_data, data[,g_var] )
		colnames( srm_data )[ncol(srm_data)] <- g_var
	}

	#- are there any predictors?
	mu_preds <- list()
	for ( i in seq_along( formulas ) ) {
		tmpX 	 <- model.matrix.lm( formulas[[i]], data, na.action = "na.pass" )
		tmpX 	 <- as.data.frame( tmpX )
		mu_preds[[i]] <- names( tmpX )
		new_cols <- setdiff( names(tmpX), colnames(srm_data) )
	    if ( length(new_cols) > 0 ) {
	        srm_data <- cbind( srm_data, tmpX[, new_cols, drop=FALSE] )
	    }
	}
	
	#- add a dyad-level identifier?
	d_var_type <- NULL
	if ( is.null( d_var ) ) {
	  	srm_data   <- srm_make_dyad_number( srm_data = srm_data, p_var = p_var, 
	  		g_var = g_var )  	
	  	d_var      <- "Dyad"
	  	d_var_type <- "Dyad_type" 
	}

	#- now everything is in wide format; we change this to long format if necessary:
	final_data_frame <- NULL
	no_var <- length( names_y )
	
	if ( no_var > 1 ) {
		
		for ( i in 1:no_var ) {
			
			tmp_data <- srm_data[,c( g_var,p_var,d_var, d_var_type, 
								     mu_preds[[i]], 
				                     names_y[i] )]
			colnames( tmp_data )[ ncol( tmp_data ) ] <- "y"
			tmp_data$measure <- i
			
			#- delete missings for relevant variables:
    		pred_cols <- mu_preds[[i]]  # inkl. "(Intercept)"
			tmp_data  <- tmp_data[ complete.cases( tmp_data[, c(pred_cols, "y")] ), ]
			final_data_frame <- rbind( final_data_frame, tmp_data )
		
		}
		
		names_y <- "y"

	} else {
		pred_cols        <- mu_preds[[1]]
		final_data_frame <- srm_data[ complete.cases( srm_data[, c(pred_cols, names_y)] ), ]
		final_data_frame$measure <- 1		
	
	}

	#- output:
	names_list <- list( g_var = g_var, p_var = p_var, d_var = d_var, d_var_type = d_var_type,
		outcome = names_y, mu_preds = mu_preds, no_var = no_var )
	out <- list( srm_data = final_data_frame, names_list = names_list )
	return( out )

}

