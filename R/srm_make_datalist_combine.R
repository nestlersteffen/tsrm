
#---- this function combines the design matrices; prior to that columns are added to group-specific matrices:

srm_make_datalist_combine <- function( data_list_groups = NULL, no_var = NULL )
{

	#- get infos:
	tmp_list  <- data_list_groups[["data_list"]]
	groupinfo <- data_list_groups[["groupinfo"]]
	ngroups   <- nrow( groupinfo )

	#- information to build big matrices
	max_groupsize <- max( groupinfo[,1] )
	max_pcolsize  <- no_var*(max_groupsize*2)#no_var*(max_groupsize*(max_groupsize+1))/2
	max_dcolsize  <- no_var*(max_groupsize*(max_groupsize-1))#2*max_pcolsize

	#- iterate across groups to add columns if necessary:
	for ( ng in seq( ngroups ) ) {

		#- when current group has size smaller than max size we adapt the matrices:
		if ( groupinfo[ng,1] < max_groupsize ) {

			#- get matrices:
			Zp <- tmp_list[[ng]][["Zp"]]
			Zd <- tmp_list[[ng]][["Zd"]]
			# Ze <- ng_data_list[[ng]][["Ze"]]
			
			#- make them greater
			Zp <- cbind( Zp, matrix( 0, nrow = nrow( Zp ), ncol = max_pcolsize - ncol( Zp ) ) )
			Zd <- cbind( Zd, matrix( 0, nrow = nrow( Zd ), ncol = max_dcolsize - ncol( Zd ) ) )

			#- now add them:
			tmp_list[[ng]][["Zp"]] <- Zp
			tmp_list[[ng]][["Zd"]] <- Zd
			#tmp_list[[ng]][["Ze"]] <- Ze

		}

	}

	#- combine all matrices into big ones:
	y  <- do.call( "rbind", lapply( tmp_list, function(x) x$y ) )
	X  <- do.call( "rbind", lapply( tmp_list, function(x) x$X ) )
	Zg <- do.call( "rbind", lapply( tmp_list, function(x) x$Zg ) )
	Zp <- do.call( "rbind", lapply( tmp_list, function(x) x$Zp ) )	
	Zd <- do.call( "rbind", lapply( tmp_list, function(x) x$Zd ) )	
	#Ze <- do.call( "rbind", lapply( data_list, function(x) x$Ze ) )
	np <- max( groupinfo[,1] )
	nd <- max( groupinfo[,2] )
	nv <- no_var

	#- output:
	data_list <- list( y = y, X = X, Zg = Zg, Zp = Zp, Zd = Zd, #Ze = Ze, 
		groupinfo = groupinfo, np = np, nd = nd, nv = nv )
	return( data_list )

}