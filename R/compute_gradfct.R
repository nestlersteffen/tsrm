
#---- this function computes the log-likelihood across all groups: 

compute_gradfct <- function( parm=NULL, parm_list=NULL, parm_table=NULL,
	data_list=NULL, args_list=NULL, model=c("srm","tsrm"), both=TRUE ) 
{
	#- get model-specific functions.
	model <- match.arg( model )
	if ( model == "srm" ) {
		include_free_parms <- srm_include_free_parameters
        sigma_derivatives  <- srm_sigma_derivatives
	} else {
		include_free_parms <- tsrm_include_free_parameters
        sigma_derivatives  <- tsrm_sigma_derivatives
	}

	#- insert parm
	parm_list <- include_free_parms( parm=parm, parm_list=parm_list, 
		parm_table=parm_table )

	#- get covariance matrices
	parm_list <- get_sigmas( parm_list=parm_list, model=model )

	#- get parm matrices:
	BETA 	<- parm_list[["BETA"]]
	SIGMA_G <- parm_list[["SIGMA_G"]]
	SIGMA_P <- parm_list[["SIGMA_P"]]
	SIGMA_D <- parm_list[["SIGMA_D"]]
    SIGMA_T <- parm_list[["SIGMA_T"]]

	#- get matrices and groupinfo:
	groupinfo <- data_list[["groupinfo"]]
	y  <- data_list[["y"]]
	X  <- data_list[["X"]]
	Zp <- data_list[["Zp"]]
	Zd <- data_list[["Zd"]]
	Zg <- data_list[["Zg"]] 
	Zt <- data_list[["Zt"]] 
	np <- data_list[["np"]]
	nd <- data_list[["nd"]]
	nt <- data_list[["nt"]]
	nv <- data_list[["nv"]]
	Pp <- data_list[["Pp"]]
	Pd <- data_list[["Pd"]]
	Pt <- data_list[["Pt"]]
	
	#- make the matrices "big"
	SIGMA_P <- Pp %*% ( diag(1,np) %x% SIGMA_P ) %*% t( Pp )
	SIGMA_D <- Pd %*% ( diag(1,nd) %x% SIGMA_D ) %*% t( Pd )
	SIGMA_T <- if ( ncol(Zt) > 0 ) Pt %*% ( diag(1, nt) %x% parm_list[["SIGMA_T"]] ) %*% t( Pt ) else SIGMA_T

	#- get no. groups:
	ngroups  <- nrow( groupinfo )
	gradient <- matrix( 0, nrow = ngroups, ncol = max( parm_table$index) ) 
	llfct    <- rep( 0, ngroups )

	#- let's go:
	for ( ng in seq( ngroups ) ) {
			
		#- get indices:
		idx1 <- groupinfo[ng,5] - groupinfo[ng,4] + 1 
		idx2 <- groupinfo[ng,5]

		#- fixed groups?
		BETA_ng <- BETA
		parm_table_ng <- parm_table
		if ( args_list$fixed_group ) {
			if ( length( BETA ) > ngroups ) {
				idx <- c( ng, (ngroups+1):length(BETA) )
			} else {
				idx <- ng 
			}
			BETA_ng <- as.matrix( BETA[idx,1] )
			if ( length( which( parm_table_ng$type == "BETA" ) ) > 0 ) {
				parm_table_ng <- parm_table[ ( parm_table$type == "BETA" & parm_table$pos1 %in% idx ) | parm_table$type != "BETA", ]
				parm_table_ng[ parm_table_ng$type == "BETA", ]$pos1 <- seq( 1, nrow( BETA_ng ), 1 )
			}	
		}
			
		#- get group and triadic matrices:
		Zg_ng <- if (ncol(Zg) > 0) as.matrix( Zg[idx1:idx2, ] ) else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)
		Zt_ng <- if (ncol(Zt) > 0) Zt[idx1:idx2, ] else matrix(0, nrow = idx2 - idx1 + 1, ncol = 0)

		if ( args_list$use_rcpp && !args_list$large ) {

			#- get position matrix and type in case of rcpp:
			parm_mat  <- as.matrix( parm_table_ng[,c("pos1","pos2","ntype","index")] )
			storage.mode( parm_mat ) <- "integer"
			
			#- compute ...
			tmp_grad <- gradient_singlegroup_export( parm_list=parm_list, 
				ty=as.matrix(y[idx1:idx2]), tX=as.matrix(X[idx1:idx2,]), 
				tZg=Zg_ng, tZp=Zp[idx1:idx2,], tZd=Zd[idx1:idx2,], tZt=Zt_ng,
				Pp=Pp, Pd=Pd, Pt=Pt,
				SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, SIGMA_T=SIGMA_T,
				BETA=as.vector(BETA_ng), np=np, nd=nd, nt=nt, parm_mat=parm_mat, 
				with_reml=args_list$with_reml, model=model )

		} else if ( args_list$use_rcpp && args_list$large ) {

			#- make sparse matrices:
			tZp     <- as(Zp[idx1:idx2,], "dgCMatrix")
			tZd     <- as(Zd[idx1:idx2,], "dgCMatrix")
			Zt_ng   <- as(Zt_ng, "dgCMatrix")
			Pp      <- as(Pp, "dgCMatrix")
			Pd      <- as(Pd, "dgCMatrix")
			Pt      <- as(Pt, "dgCMatrix")
			SIGMA_G <- if (model == "tsrm") matrix(0, 0, 0) else SIGMA_G
			SIGMA_P <- as(SIGMA_P, "dgCMatrix")
			SIGMA_D <- as(SIGMA_D, "dgCMatrix")
			SIGMA_T <- if (model == "srm") as(matrix(0, 0, 0), "dgCMatrix") else {
    			SIGMA_T |> as("dgCMatrix")
			}

			#- get position matrix and type in case of rcpp:
			parm_mat  <- as.matrix( parm_table_ng[,c("pos1","pos2","ntype","index")] )
			storage.mode( parm_mat ) <- "integer"
			
			#- compute ...
			tmp_grad <- gradient_singlegroup_sparse_export( parm_list=parm_list, 
				ty=as.matrix(y[idx1:idx2]), tX=as.matrix(X[idx1:idx2,]), 
				tZg=Zg_ng, tZp=tZp, tZd=tZd, tZt=Zt_ng, Pp=Pp, Pd=Pd, Pt=Pt,
				SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, SIGMA_T=SIGMA_T,
				BETA=as.vector(BETA_ng), np=np, nd=nd, nt=nt, parm_mat=parm_mat, 
				with_reml=args_list$with_reml, model=model )

		} else {

			tmp_grad <- compute_gradient_singlegroup(
				parm_list=parm_list, parm_table=parm_table_ng, 
				np=np, nd=nd, nt=nt, 
				y=y[idx1:idx2], X=X[idx1:idx2,],
				Zg=Zg_ng, Zp=Zp[idx1:idx2, ], Zd=Zd[idx1:idx2, ], Zt=Zt_ng, Pp=Pp, Pd=Pd, Pt=Pt,
	            BETA=BETA_ng, SIGMA_G=SIGMA_G, SIGMA_P=SIGMA_P, SIGMA_D=SIGMA_D, SIGMA_T=SIGMA_T, 
				args_list=args_list, sigma_derivatives=sigma_derivatives )

		}

		llfct[ng]     <- tmp_grad[ 1]
		gradient[ng,] <- tmp_grad[-1]

	}
	
	#- output:
	if ( both ) {
		return( list( "objective"=-1*sum(llfct), "gradient"=-1*colSums(gradient) ) )
	} else {
		return( -1*colSums(gradient) )
	}

}