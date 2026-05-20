
#---- this function generates the design matrices per group for the tsrm:

tsrm_make_dataframe <- function( names_y=NULL, names_X=NULL, p_var=NULL,
  	g_var=NULL, formulas=NULL, data=NULL ) 
{

	data_names <- colnames( data )

	#- some checks...
	if ( !all( p_var %in% data_names ) | is.null( p_var ) ) {
		stop("Something is wrong with the person-level identifiers.")
	}
	if ( !( all( names_y %in% data_names ) ) ) {
		stop("Something is wrong with the outcome variable.")
	}

	#- we build a dataframe:
	tmp_data <- data[,c(p_var,names_y)]
	
	#- add a group-variable if necessary:
	if ( is.null( g_var ) ) {
	  	tmp_data$group <- 1
	  	g_var <- "group"
	} else {
		tmp_data <- cbind( tmp_data, data[,g_var] )
		colnames( tmp_data )[ncol(tmp_data)] <- g_var
	}

	#- are there any predictors?
	mu_preds <- list()
	for ( i in seq_along( formulas ) ) {
		tmpX 	 <- model.matrix.lm( formulas[[i]], data, na.action = "na.pass" )
		tmpX 	 <- as.data.frame( tmpX )
		mu_preds[[i]] <- names( tmpX )
		new_cols <- setdiff( names(tmpX), colnames( tmp_data ) )
	    if ( length(new_cols) > 0 ) {
	        tmp_data <- cbind( tmp_data, tmpX[, new_cols, drop=FALSE] )
	    }
	}
	
	#- add a dyad-level identifier:
	tmp <- c("a","b","c")
	for ( i in 1:2 ) {
		for ( j in (i+1):3 ) {
			tmp_data <- make_dyad_number( data=tmp_data, p_var=p_var[c(i,j)], 
	        	g_var=g_var, dyad_name=paste0( tmp[i],tmp[j]), maxg = 1e2 )
		}
	}
	d_var      <- c("ab","ac","bc")
	d_var_type <- paste0( d_var, "_type" )
	
	#- add a triad-level identifier:
	tmp_data   <- make_triad_number( data=tmp_data, p_var=p_var, g_var=g_var )
	t_var      <- "Triad"
	t_var_type <- "Triad_type"
	
	#- now everything is in wide format; we change this to long format if necessary:
	final_data_frame <- NULL
	no_var <- length( names_y )
	
	if ( no_var > 1 ) {
		
		for ( i in 1:no_var ) {
			
			tmp_data <- tmp_data[,c( g_var,p_var,d_var,d_var_type,
									 t_var,t_var_type,
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
		final_data_frame <- tmp_data[ complete.cases( tmp_data[, c(pred_cols, names_y)] ), ]
		final_data_frame$measure <- 1		
	
	}

	#- output:
	names_list <- list( g_var=g_var, p_var=p_var, d_var=d_var, d_var_type=d_var_type,
		t_var=t_var, t_var_type=t_var_type, outcome=names_y, mu_preds=mu_preds, 
		no_var=no_var )
	out <- list( data=final_data_frame, names_list=names_list )
	return( out )

}