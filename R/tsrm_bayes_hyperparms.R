
#---- function to make a list with hyperparameters:

tsrm_bayes_hyperparms <- function( args_list=NULL, data_list=NULL, 
	person=NULL, dyads=NULL ) 
{
	
	#- get specifications:
	mcmc_args <- args_list$mcmc_args
	y.scale   <- mcmc_args$y.scale
	mb 	      <- mcmc_args$mb
	invMb     <- mcmc_args$invMb
	IW.nup	  <- mcmc_args$IW.nup
	IW.Sp     <- mcmc_args$IW.Sp
	LKJ.eta   <- mcmc_args$LKJ.eta
	Ht.scale  <- mcmc_args$Ht.scale
	Ht.nu     <- mcmc_args$Ht.nu
	tau       <- mcmc_args$tau
	tau_sd    <- mcmc_args$tau_sd
	tau_rho   <- mcmc_args$tau_rho
	IW.nud	  <- mcmc_args$IW.nud
	IW.Sd     <- mcmc_args$IW.Sd
	IW.nut	  <- mcmc_args$IW.nut
	IW.St     <- mcmc_args$IW.St
	IG.at     <- mcmc_args$IG.at
	IG.bt     <- mcmc_args$IG.bt

	#- get scale of outcome:
	ys <- if ( is.null( y.scale ) ) sd( data_list$y ) else y.scale
		
	#- now define hyperparms:
	if ( is.null( mb ) ) {
		mcmc_args$mb <- matrix( 0, nrow = ncol( data_list$X ), ncol = 1)
	} else {
		if ( !is.matrix( mb ) ) {
			stop( "mb must a be a k x 1 - matrix (k = no. of preds)")
		}
		if ( nrow( mb ) != ncol( data_list$X ) ) {
			stop( "mb does not match with no. of preds")
		}
	}
	
	if ( is.null( invMb ) ) {
		mcmc_args$invMb <- diag( 1/( 100*ys^2 ), ncol( data_list$X ) )
	} else { 
		if ( !is.matrix( invMb ) ) {
			stop( "invMb must a be a k x k - matrix (k = no. of preds)")
		}
		if ( length( diag( invMb ) ) != ncol( data_list$X ) ) {
			stop( "invMb does not match with no. of preds")
		}
	}

	#%%%%%%%%% SIGMA_P

	if ( is.null( IW.nup ) ) IW.nup <- 5
	mcmc_args$IW.nup <- IW.nup

	if ( is.null( IW.Sp ) ) {
		mcmc_args$IW.Sp <- (IW.nup-4)*diag( ys^2 / 4, 3 )
	} else {
		if ( !is.matrix( IW.Sp ) ) {
			stop( "IW.Sp must a be a 3 x 3 - matrix")
		}
		if ( !all( dim( IW.Sp ) == c( 3, 3 ) ) ) {
			stop( "IW.Sp is not a 3 x 3 matrix")
		}
	}

	if ( is.null( mcmc_args$LKJ.eta ) ) {
		mcmc_args$LKJ.eta <- 1
	}

	if ( is.null( mcmc_args$Ht.scale ) ) {
		mcmc_args$Ht.scale <- ys
	}

	if ( is.null( mcmc_args$Ht.nu ) ) {
		mcmc_args$Ht.nu <- 3
	}

	if ( is.null( mcmc_args$tau ) ) {
		mcmc_args$tau <- 0.05
	}

	if ( is.null( mcmc_args$tau_sd ) ) {
		mcmc_args$tau_sd <- 0.20
	}

	if ( is.null( mcmc_args$tau_rho ) ) {
		mcmc_args$tau_rho <- 0.50
	}

	#%%%%%%% SIGMA_D

	if ( is.null( IW.nud ) ) IW.nud <- 5
	mcmc_args$IW.nud <- IW.nud

	if ( is.null( IW.Sd ) ) {
		mcmc_args$IW.Sd <- list( Sd1=(IW.nud-4)*diag( ys^2 / 4, 3), Sd2=(IW.nud-4)*diag( ys^2 / 4, 3 ) )
	} else {
		if ( !is.list( IW.Sd ) ) {
			stop( "IW.Sd must be a list of two 3 x 3 matrices.")
		} else {
			if ( !is.matrix( IW.Sd[[1]] ) | !is.matrix( IW.Sd[[2]] ) ) {
				stop( "One element in IW.Sd is not a 3 x 3 - matrix")
			}
			if ( !all( dim( IW.Sd[[1]] ) == c( 3, 3 ) ) || !all( dim( IW.Sd[[2]] ) == c( 3, 3 ) ) ) {
				stop( "One element in IW.Sd is not a 3 x 3 matrix")
			}
		}	
	}

	#%%%%%%% SIGMA_T
	
	if ( is.null( IW.nut ) ) IW.nut <- 4
	mcmc_args$IW.nut <- IW.nut

	if ( is.null( IW.St ) ) {
		mcmc_args$IW.St <- (IW.nut-3)*diag( ys^2 / 4, 2 )
	} else {
		if ( !is.matrix( IW.St ) ) {
			stop( "Per construction IW.St must a be a 2 x 2 - matrix")
		}
		if ( !all( dim( IW.St ) == c( 2, 2 ) ) ) {
			stop( "IW.St is not a 2 x 2 matrix")
		}
	}

	if ( is.null( mcmc_args$IG.at ) ) {
		mcmc_args$IG.at <- 1
	}

	if ( is.null( mcmc_args$IG.bt ) ) {
		mcmc_args$IG.bt <- 1
	}

	#- return
	args_list$mcmc_args <- mcmc_args
	return( args_list )
}