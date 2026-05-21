
make_args_list <- function( 
	use_rcpp = TRUE,
	with_reml = TRUE,
	random_group = FALSE,
	fixed_group = FALSE,
	type_ses = "Standard",
	method = "c++",
	n_threads = 1L,
	large = FALSE,
	maxit = 1000L,
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
  
	if ( !( method %in% c("c++","tmb","Bayes") ) ) {
		stop("Undefined method defined.")
	}

	if( !is.list( mcmc_args ) ) { 
		stop("'mcmc_args' has to be a list.") 
	}

	if( !is.list( nlminb_list ) ) {
		stop("'nlminb_list' has to be a list.") 
	}

	if ( fixed_group ) {
		random_group <- FALSE
	}

	if ( random_group ) {
		fixed_group <- FALSE
	}

	#- make args for ML estimation:
	if ( method == "c++" | method == "tmb" ) {
		#- make default list:
		default_list <- list( 
			iter.max = maxit,
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
	if ( method == "Bayes" ) {
		#- make default list:
		default_list <- list( mb = NULL, invMb = NULL, 
			# all things for sigmap
			prior_sigmap = "IW", # or LKJ
			IW.nu0 = NULL, IW.S0 = NULL,
			LKJ.eta = NULL, LKJ.scale = NULL, LKJ.nu = NULL, LKJ.tau = NULL,
			# all things for sigmad
			sigmad_update = "joint", # vs. seperate
			prior_sigmad = "HalfCauchy", # vs. Halft vs. Uniform
			HCd_s = NULL, Htd_scale = NULL, Htd_nu = NULL,
			tau_sd_d = NULL, tau_rho_d = NULL,
			# all thing for sigmat:
			sigmat_update = "joint", # vs. separate
			prior_sigmat = "HalfCauchy", # vs. Halft vs. Uniform
			HCt_s = NULL, Htt_scale = NULL, Htt_nu = NULL,
			tau_sd_t = NULL, tau_rho_t = NULL,
			verbose = verbose )
		#- add to final list:
		mcmc_args <- c( mcmc_args, 
			default_list[ !( names( default_list ) %in% names( mcmc_args ) ) ] )
	}
	
	args_list <- list( 
		use_rcpp = use_rcpp, # use rcpp?
		with_reml = with_reml,
		random_group = random_group,
		fixed_group = fixed_group,
		method = method,
		n_threads = 1L,
		large = large,
		type_ses = type_ses,
		nlminb_ctrl = nlminb_ctrl,
		mcmc_args = mcmc_args,
		raneffs_type = raneffs_type,
		no_pvs = no_pvs, )
	return( args_list )
} 