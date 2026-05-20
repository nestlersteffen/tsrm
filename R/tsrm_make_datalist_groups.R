
#---- this function generates the design matrices per group:

tsrm_make_datalist_groups <- function( data=NULL, names_list=NULL, args_list=NULL ) 
{

	#- get names:
	g_var    <- names_list[["g_var"]]
	p_var    <- names_list[["p_var"]]
	d_var    <- names_list[["d_var"]]
	t_var    <- names_list[["t_var"]]
	out      <- names_list[["outcome"]]
	mu_preds <- names_list[["mu_preds"]] 

	#- rr-group information:
	groups    <- unique( data[,g_var] )
	ngroups   <- base::length( groups)
	tmp_list  <- vector( "list", ngroups )
	groupinfo <- matrix( 0, nrow = ngroups, ncol = 5 )

	#- how many variables?
	no_var    <- names_list[["no_var"]] 

	#- Step 1: iterate through the groups and each measure to build the matrices:
	for ( ng in seq(ngroups) ) {

		tmp_y <- tmp_X <- tmp_Zg <- tmp_Zp <- tmp_Zd <- tmp_Zt <- vector("list",no_var)

		for ( nv in seq( no_var ) ) {

			#- temporary data frame for rr-variables for the measure 
			tmp_data <- data[data[,g_var] == groups[ng] & data[,"measure",] == nv,]

			#- get outcome vector and predictors:
			tmp_y[[nv]] <- as.matrix( tmp_data[,out] )
			tmp_X[[nv]] <- as.matrix( tmp_data[,mu_preds[[nv]]] )

			#- ---- make group design matrices (just a zero-column matrix )

			tmp_Zg[[nv]] <- matrix( 0, nrow=nrow( tmp_data ), ncol=0 )

			#- ---- make person design matrix

			#- get the person identifiers:
			persons    <- sort( unique( c( tmp_data[,p_var[1]], tmp_data[,p_var[2]], tmp_data[,p_var[3]] ) ) )
			no_persons <- length( persons )
		  
		  	#- build matrix:
			Zp <- matrix( 0, ncol = 3*no_persons, nrow = nrow( tmp_data ) )
			for ( ii in 1:no_persons ) {
				idx_actor   <- which( tmp_data[,p_var[1]]==persons[ii] )
				idx_partner <- which( tmp_data[,p_var[2]]==persons[ii] )
				idx_judge   <- which( tmp_data[,p_var[3]]==persons[ii] )
				Zp[idx_actor  ,3*ii-2] <- 1
			    Zp[idx_partner,3*ii-1] <- 1
			    Zp[idx_judge  ,3*ii]   <- 1
			}
			tmp_Zp[[nv]] <- Zp
			  
			#- ------ make dyad effect design matrix

			dyads    <- sort( unique( c( tmp_data[,d_var[1]], tmp_data[,d_var[2]], tmp_data[,d_var[3]] ) ) ) 
			no_dyads <- length( unique( c( tmp_data[,d_var[1]], tmp_data[,d_var[2]], tmp_data[,d_var[3]] ) ) )

			Zd <- matrix(0, ncol = 6*no_dyads, nrow = nrow(tmp_data))
			for (ii in 1:nrow(tmp_data)) {
			    # ab-Effekt
			    dpos_ab <- which(dyads == tmp_data[ii, d_var[1]])
			    offset_ab <- ifelse(tmp_data[ii, "ab_type"] == 1, 0, 1)
			    Zd[ii, 6*(dpos_ab - 1) + 1 + offset_ab] <- 1
			    # ac-Effekt
			    dpos_ac <- which(dyads == tmp_data[ii, d_var[2]])
			    offset_ac <- ifelse(tmp_data[ii, "ac_type"] == 1, 0, 1)
			    Zd[ii, 6*(dpos_ac - 1) + 3 + offset_ac] <- 1
			    # bc-Effekt
			    dpos_bc <- which(dyads == tmp_data[ii, d_var[3]])
			    offset_bc <- ifelse(tmp_data[ii, "bc_type"] == 1, 0, 1)
			    Zd[ii, 6*(dpos_bc - 1) + 5 + offset_bc] <- 1
			}
			tmp_Zd[[nv]] <- Zd

			#- ------ make triad design matrix
			triads <- sort( unique( tmp_data[,"Triad"] ) ) 
			no_triads <- length( triads )

			Zt <- matrix( 0, ncol = 6*no_triads, nrow = nrow( tmp_data ) )
			for ( ii in 1:nrow( tmp_data ) ) { 
			    # triad number and type of triad: 
			    trval <- tmp_data[ii, "Triad"]
			    trtyp <- tmp_data[ii, "Triad_type"]
			    # triad is the xx triad in triads:
			    trpos <- which( triads == trval )
			    # compute colindex:
			    col_index <- (trpos - 1) * 6 + trtyp
			    Zt[ii, col_index] <- 1
			}
			tmp_Zt[[nv]] <- Zt

		} 
		
		#- make final matrices:
		y  <- do.call( "rbind", tmp_y )
		X  <- make_bldiag( tmp_X )
		Zg <- make_bldiag( tmp_Zg )
		Zp <- make_bldiag( tmp_Zp )
		Zd <- make_bldiag( tmp_Zd )
		Zt <- make_bldiag( tmp_Zt )

		#- ------ save everything in the list
		tmp_list[[ng]] <- list( y=y, X=X, Zg=Zg, Zp=Zp, Zd=Zd, Zt=Zt )

		#- ------ save the rest in the matrix:
		groupinfo[ng,1] <- no_persons
		groupinfo[ng,2] <- no_dyads
		groupinfo[ng,3] <- no_triads
		groupinfo[ng,4] <- base::length( y )
	
	} # end for ng

	groupinfo[,5] <- base::cumsum( groupinfo[,4] )

	#- output:
	out <- list( data_list=tmp_list, groupinfo=groupinfo ) 
  	return( out )

}