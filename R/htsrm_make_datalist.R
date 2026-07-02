
#---- this function generates the design matrices per group:

htsrm_make_datalist <- function( data=NULL, names_list=NULL, args_list=NULL ) 
{

	#- STEP 1: we have a dataframe with two variables, measure = 1 is the triadic round-robin variable and 
	#  measure = 2 is the dyadic variable; for each measure, we obtain the design-matrices

	#- start with measure 1, the triadic variable
	tsrm_data      <- data[data[,"measure",] == 1,]
	tmp_names_list <- names_list
	tmp_names_list$no_var <- 1
	tsrm_datalist_groups  <- tsrm_make_datalist_groups( data=tsrm_data, names_list=tmp_names_list, args_list=args_list )

	#- start with measure 2, the dyadic variable:
	srm_data             <- data[data[,"measure",] == 2,]
	srm_data$measure     <- 1
	srm_data$Dyad_type   <- srm_data$ab_type
	tmp_names_list$p_var <- tmp_names_list$p_var[1:2]
	tmp_names_list$d_var <- "ab"
	tmp_names_list$mu_preds[[1]] <- tmp_names_list$mu_preds[[2]]
	srm_datalist_groups  <- srm_make_datalist_groups( data=srm_data, names_list=tmp_names_list, args_list=args_list  )

	#- STEP 2: we pad the design matrices within the groups to max_groupsize
	tsrm_datalist_groups <- htsrm_pad_to_max( data_list_groups=tsrm_datalist_groups, model="tsrm" )
	tsrm_datalist  <- tsrm_datalist_groups$data_list
	tsrm_groupinfo <- tsrm_datalist_groups$groupinfo

	srm_datalist_groups  <- htsrm_pad_to_max( data_list_groups=srm_datalist_groups, model="srm" )
	srm_datalist  <- srm_datalist_groups$data_list
	srm_groupinfo <- srm_datalist_groups$groupinfo

	#- STEP 3: we iterate through the groups and merge the design matrices:
	g_var     <- names_list[["g_var"]]
	groups    <- unique( data[,g_var] )
	ngroups   <- base::length( groups)
	tmp_list  <- vector( "list", ngroups )
	groupinfo <- matrix( 0, nrow = ngroups, ncol = 5 )

	for ( ng in seq(ngroups) ) {

		#- get the group specific design matrices for the two variables
		ng_tsrm <- tsrm_datalist[[ng]]
		ng_srm  <- srm_datalist[[ng]]

		#- how many triadic and dyadic judgments:
		n1 <- nrow( ng_tsrm$y ); n2 <- nrow( ng_srm$y )

		#- now build the matrices:
		y  <- rbind( ng_tsrm$y, ng_srm$y )
		X  <- make_bldiag( list( ng_tsrm$X,  ng_srm$X  ) )
        Zg <- matrix( 0, nrow = n1 + n2, ncol = 0 )
        Zp <- make_bldiag( list( ng_tsrm$Zp, ng_srm$Zp ) )   # 3*np1  | 2*np2
        Zd <- make_bldiag( list( ng_tsrm$Zd, ng_srm$Zd ) )   # 6*nd1  | 2*nd2
		Zt <- rbind( ng_tsrm$Zt, matrix( 0, nrow = n2, ncol = ncol( ng_tsrm$Zt ) ) )

		#- ------ save everything in the list
		tmp_list[[ng]] <- list( y=y, X=X, Zg=Zg, Zp=Zp, Zd=Zd, Zt=Zt )

		#- ------ save the rest in the matrix: # only balanced Designs?
		groupinfo[ng,1] <- max( c( tsrm_groupinfo[ng,1], srm_groupinfo[ng,1] ) )  # np
		groupinfo[ng,2] <- max( c( tsrm_groupinfo[ng,2], srm_groupinfo[ng,2] ) )  # nd
		groupinfo[ng,3] <- tsrm_groupinfo[ng,3] 
		groupinfo[ng,4] <- base::length( y )
	
	} # end for ng

	groupinfo[,5] <- base::cumsum( groupinfo[,4] )

	#- STEP 4: we combine alle group specific matrices:
	y  <- do.call( "rbind", lapply( tmp_list, function(x) x$y ) )
	X  <- do.call( "rbind", lapply( tmp_list, function(x) x$X ) )
	Zg <- do.call( "rbind", lapply( tmp_list, function(x) x$Zg ) )
	Zp <- do.call( "rbind", lapply( tmp_list, function(x) x$Zp ) )	
	Zd <- do.call( "rbind", lapply( tmp_list, function(x) x$Zd ) )	
	Zt <- do.call( "rbind", lapply( tmp_list, function(x) x$Zt ) )
	np <- max( groupinfo[,1] )
	nd <- max( groupinfo[,2] )
	nt <- max( groupinfo[,3] )
	nv <- 2

	#- STEP 5: get the permutation matrix:
	tmp_p <- make_permutationmatrix( n_units=np, block_size=c(3,2), nv=nv )
	tmp_d <- make_permutationmatrix( n_units=nd, block_size=c(6,2), nv=nv )
	tmp_t <- make_permutationmatrix( n_units=nt, block_size=c(6,0), nv=nv )

	#- output:
	data_list <- list( y=y, X=X, Zg=Zg, Zp=Zp, Zd=Zd, Zt=Zt, groupinfo=groupinfo, 
		np=np, nd=nd, nt=nt, nv=nv, Pp=tmp_p$P, Pd=tmp_d$P, Pt=tmp_t$P,
		perm_p=tmp_p$perm, perm_d=tmp_d$perm, perm_t=tmp_t$perm )
	return( data_list )

}