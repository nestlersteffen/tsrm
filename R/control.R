
make_args_list <- function( 
	use_rcpp = TRUE,
	method = "reml", # ml, cml, creml, bayes 
	with_reml = TRUE,
	random_group = FALSE,
	fixed_group = FALSE,
	constrained = FALSE,
	type_ses = "Standard",
	n_threads = 1L,
	starts = "none",
	large = FALSE,
	maxit = 1000L,
	maxeval = 5000L,
	abs_tol = (.Machine$double.eps*10),
	rel_tol = 1e-9,
	step_min = 2.2e-10,
	x_tol = 1.5e-8,
	verbose = TRUE,
	nlminb_list = list(),
	mcmc_args = list(),
	raneffs_type = "EB",
	no_pvs = 50,
	...) 
{
  
	if ( !( method %in% c("ml","reml", "cml", "creml", "bayes") ) ) {
		stop("Undefined value for 'method' used.")
	}

	if ( !( starts %in% c("none","anova") ) ) {
		stop("Undefined value for 'starts' used.")
	}

	if( !is.list( mcmc_args ) ) { 
		stop("'mcmc_args' has to be a list.") 
	}

	if( !is.list( nlminb_list ) ) {
		stop("'nlminb_list' has to be a list.") 
	}

	if ( method %in% c( "ml", "cml" ) ) with_reml <- FALSE
	if ( method %in% c( "cml", "creml" ) ) constrained <- TRUE

	if ( fixed_group ) random_group <- FALSE
	if ( random_group ) fixed_group <- FALSE

	# if ( !large | !use_rcpp ) {
	# 	use_rcpp <- FALSE
	# 	large    <- FALSE
	# }

	#- make args for ML estimation:
	if ( method %in% c("ml", "reml", "cml", "creml" ) ) {
		#- make default list:
		default_list <- list( 
			iter.max = maxit,
			eval.max = maxeval,
			abs.tol = abs_tol, 
			rel.tol = rel_tol,
			step.min = step_min, 
			x.tol = x_tol,
			trace = 0 
		)
		#- check verbose:
		if ( verbose ) { default_list$trace = 1 }
		
		#- add to final list:
		nlminb_ctrl <- c( nlminb_list, 
        	default_list[ !( names( default_list ) %in% names( nlminb_list ) ) ] )
	}

	#- make args for Bayes:
	if ( method == "bayes" ) {
	 	#- make default list:
		default_list <- list( 
			# scale of the outome:
			y.scale = NULL,
			# all things for Beta
			mb = NULL, invMb = NULL, 
			# all things for sigmap
			prior_sigmap = "IW", # or LKJ or LKJ_single
			IW.nup = NULL, IW.Sp = NULL,
			LKJ.eta = NULL, Ht.scale = NULL, Ht.nu = NULL, tau = NULL,
			tau_sd = NULL, tau_rho = NULL,
			# all things for sigmad
			IW.nud = NULL, IW.Sd = NULL,
			# all thing for sigmat:
			IW.nut = NULL, IW.St = NULL,
			IG.at = NULL, IG.bt = NULL,
			# all things for generating chain-specific start values
			sd_log=0.5, eta=1, nu=8 )
		#- add to final list:
		mcmc_args <- c( mcmc_args, 
			default_list[ !( names( default_list ) %in% names( mcmc_args ) ) ] )
	}
	
	args_list <- list( 
		use_rcpp = use_rcpp, # use rcpp?
		method = method,
		with_reml = with_reml,
		random_group = random_group,
		fixed_group = fixed_group,
		constrained = constrained,
		n_threads = 1L,
		starts = starts,
		large = large,
		type_ses = type_ses,
		nlminb_ctrl = nlminb_ctrl,
		mcmc_args = mcmc_args,
		raneffs_type = raneffs_type,
		no_pvs = no_pvs )
	return( args_list )
} 