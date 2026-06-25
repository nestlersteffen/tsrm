
#---- this function combines the design matrices; prior to that columns are 
#     added to group-specific matrices

htsrm_pad_to_max <- function( data_list_groups=NULL, no_var=1, model=NULL )
{

	#- get infos:
	tmp_list  <- data_list_groups[["data_list"]]
	groupinfo <- data_list_groups[["groupinfo"]]
	ngroups   <- nrow( groupinfo )

	#- get multiplier based on model:
	pmult <- 2; dmult <- 1
	psize <- 2; dsize <- 2
	if ( model == "tsrm" ) {
		pmult <- 3; dmult <- 3
		psize <- 3; dsize <- 6
	}

	#- information to build big matrices
	max_groupsize <- max( groupinfo[,1] )
	max_pcolsize  <- no_var*(max_groupsize*pmult)#no_var*(max_groupsize*(max_groupsize+1))/2
	max_dcolsize  <- no_var*(max_groupsize*(max_groupsize-1)*dmult)#2*max_pcolsize
	max_tcolsize  <- no_var*(max_groupsize*(max_groupsize-1)*(max_groupsize-2))

	#- iterate across groups to add columns if necessary:
	for ( ng in seq( ngroups ) ) {

	 	#- when current group has size smaller than max size we adapt the matrices:
	 	if ( groupinfo[ng,1] < max_groupsize ) {

	 		#- get matrices:
	 		Zp <- tmp_list[[ng]][["Zp"]]
	 		Zd <- tmp_list[[ng]][["Zd"]]
	 		Zt <- tmp_list[[ng]][["Zt"]]
			
	 		#- make them greater
	 		Zp <- cbind( Zp, matrix( 0, nrow = nrow( Zp ), ncol = max_pcolsize - ncol( Zp ) ) )
	 		Zd <- cbind( Zd, matrix( 0, nrow = nrow( Zd ), ncol = max_dcolsize - ncol( Zd ) ) )
	 		if ( ncol( Zt ) > 0 ) {
	 			Zt <- cbind( Zt, matrix( 0, nrow = nrow( Zt ), ncol = max_tcolsize - ncol( Zt ) ) )
	 		}

	 		#- add matrices to list:
	 		tmp_list[[ng]][["Zp"]] <- Zp
	 		tmp_list[[ng]][["Zd"]] <- Zd
	 		tmp_list[[ng]][["Zt"]] <- Zt

		}

	}

	return( list( data_list=tmp_list, groupinfo=groupinfo ) )

}