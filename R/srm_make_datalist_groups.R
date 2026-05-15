
#---- this function generates the design matrices per group:

srm_make_datalist_groups <- function( srm_data = NULL, names_list = NULL, 
	args_list = NULL ) 
{

	#- get names:
	g_var    <- names_list[["g_var"]]
	p_var    <- names_list[["p_var"]]
	d_var    <- names_list[["d_var"]]
	out      <- names_list[["outcome"]]
	mu_preds <- names_list[["mu_preds"]] 

	#- rr-group information:
	groups    <- unique( srm_data[,g_var] )
	ngroups   <- base::length( groups)
	tmp_list  <- vector( "list", ngroups )
	groupinfo <- matrix( 0, nrow = ngroups, ncol = 4 )

	#- how many variables?
	no_var    <- names_list[["no_var"]] 

	#- Step 1: iterate through the groups and each measure to build the matrices:
	for ( ng in seq(ngroups) ) {

		tmp_y <- tmp_X <- tmp_Zg <- tmp_Zp <- tmp_Zd <- vector("list",no_var)

		for ( nv in seq( no_var ) ) {

			#- temporary data frame for rr-variables for the measure 
			tmp_data <- srm_data[srm_data[,g_var] == groups[ng] & srm_data[,"measure",] == nv,]

			#- sort the data:
			# tmp_data1 <- subset( tmp_data, tmp_data[,p_var[1]] < tmp_data[,p_var[2]] )
			# tmp_data1 <- tmp_data1[ order(tmp_data1[,p_var[1]],tmp_data1[,p_var[2]],tmp_data1[,d_var]),]
			# tmp_data2 <- subset( tmp_data, tmp_data[,p_var[1]] > tmp_data[,p_var[2]] )
			# tmp_data2 <- tmp_data2[ order(tmp_data2[,p_var[2]],tmp_data2[,p_var[1]],tmp_data2[,d_var]),]
			# tmp_data  <- rbind( tmp_data1,tmp_data2)
			# #print( tmp_data )
				
			#- get outcome vector and predictors:
			tmp_y[[nv]] <- as.matrix( tmp_data[,out] )
			tmp_X[[nv]] <- as.matrix( tmp_data[,mu_preds[[nv]]] )

			#- ---- make group design matrices

			tmp_Zg[[nv]] <- matrix( 1.0, ncol = 1, nrow = nrow( tmp_data ) )
			  
			#- -------

			#- ---- make person design matrices

			#- get the person identifiers:
			persons    <- sort( unique( c(tmp_data[,p_var[1]], tmp_data[,p_var[2]] ) ) )
			no_persons <- length( persons )
			  
			# #- make two design matrices:
			#   Z1 <- Z2 <- matrix( 0, ncol = no_persons, nrow = nrow( tmp_data ) )
			#   for ( ii in 1:no_persons ) {
			#       idx_actor   <- which(tmp_data[,p_var[1]]==persons[ii])
			#       idx_partner <- which(tmp_data[,p_var[2]]==persons[ii])
			#       Z1[idx_actor,ii] <- 1
			#       Z2[idx_partner,ii] <- 1
			#   }
			#   tmp_Zp[[nv]] <- cbind( Z1, Z2 ) 
			Zp <- matrix( 0, ncol = 2*no_persons, nrow = nrow( tmp_data ) )
			for ( ii in 1:no_persons ) {
			 	idx_actor   <- which( tmp_data[,p_var[1]]==persons[ii] )
			 	idx_partner <- which( tmp_data[,p_var[2]]==persons[ii] )
			 	Zp[idx_actor, 2*ii-1] <- 1.0 # actor-Effekte in die ungerade Spalte
			 	Zp[idx_partner, 2*ii] <- 1.0
			}
			tmp_Zp[[nv]] <- Zp
			  
			#- -------

			#- ------ make dyad design matrices

			#- get the dyad identifier:
			dyads    <- sort( unique( c(tmp_data[,d_var] ) ) )
			no_dyads <- length( dyads )

			#- make the design matrix:
			Zd <- matrix(0, ncol = 2*no_dyads, nrow = nrow( tmp_data ) )
			# for ( dd in 1:no_dyads ) {
			#     idx_dyad <- which( tmp_data[,d_var[1]] == dyads[ dd ] )
			#     for ( zz in 1:length( idx_dyad ) ) {
			#     	Zd[idx_dyad[zz],idx_dyad[zz]] <- 1
			#     }
			# }
			for ( ii in 1:nrow(tmp_data)) {
  				dyad_value    <- tmp_data[ii, d_var]
  				dyad_type     <- tmp_data[ii, "Dyad_type"]
   				dyad_position <- which(dyads == dyad_value)
  				if (dyad_type == 1) {
					col_index <- 2*dyad_position-1#dyad_position #  
				} else {
					col_index <- 2*dyad_position #no_dyads + dyad_position # 
				}
  				Zd[ii, col_index] <- 1.0
			}
			tmp_Zd[[nv]] <- Zd
			#print( tmp_Zd )

			#- --------

			#- ------ error matrix, in case of multiple measures
			
			#Ze <- NULL
			# if ( no_vars > 1 ) {
			# 	# ....
			# }
		
		} 
		
		#- make final matrices:
		y  <- do.call( "rbind", tmp_y )
		X  <- srm_make_bldiag( tmp_X )
		Zg <- srm_make_bldiag( tmp_Zg )
		Zp <- srm_make_bldiag( tmp_Zp )
		Zd <- srm_make_bldiag( tmp_Zd )
		#Ze <- srm_make_bldiag( tmp_Ze )

		#- ------ save everything in the list
		tmp_list[[ng]] <- list( y = y, X = X, Zg = Zg, Zp = Zp, Zd = Zd )#, Ze = Ze )

		#- ------ save the rest in the matrix:
		groupinfo[ng,1] <- no_persons
		groupinfo[ng,2] <- no_dyads
		groupinfo[ng,3] <- base::length( y )
	
	} # end for ng

	groupinfo[,4] <- base::cumsum( groupinfo[,3] )

	#- output:
	out <- list( data_list = tmp_list, groupinfo = groupinfo ) 
  	return( out )

}