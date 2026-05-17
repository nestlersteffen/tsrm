
srm_control <- function( 
  	use_rcpp = TRUE, # use rcpp?
	with_reml = TRUE,
	random_group = FALSE,
	fixed_group = FALSE,
	large = FALSE,
	type_ses = "Standard",
	method = "c++",
	verbose = FALSE,
	gr_maxit = 1000L,
	gr_abs_tol = (.Machine$double.eps*10),
	gr_rel_tol = 1e-9,
	gr_step_min = 2.2e-10,
	gr_x_tol = 1.5e-8,
	nlminb_list = NULL,
	with_ses = TRUE, 
	raneffs_type = "EB",
	no_pvs = 50,
	mcmc_args = list(),
	... ) 
{
  
	nlminb_ctrl <- NULL

	if ( !( method %in% c("c++","tmb") ) ) {
		stop("Undefined method defined.")
	}

	if( !is.list( mcmc_args ) ) { 
		stop("'mcmc_args' has to be a list.") 
	}

	if ( is.null( nlminb_list ) ) {
		nlminb_ctrl <- list( 
	    	iter.max = gr_maxit,
	    	abs.tol = gr_abs_tol, 
	    	rel.tol = gr_rel_tol,
	    	step.min = gr_step_min, 
	    	x.tol = gr_x_tol,
	    	trace = 0 )
		if ( verbose ) nlminb_ctrl$trace = 1
	}

	#if ( method == "Bayes" ) {
        #- make default list:
        default_list <- list( mb = NULL, invMb = NULL, 
   				# all things for sigmap
   				prior_sigmap = "IW", # or LKJ
   				IW.nu0 = NULL, IW.S0 = NULL, 
   				LKJ.eta = NULL, LKJ.scale = NULL, LKJ.nu = NULL, LKJ.tau = NULL,
   				# all things for sigmas
   				prior_sigmad = "HalfCauchy", # or Halft
   				HC.s = NULL, Ht.scale = NULL, Ht.nu = NULL, 
   				tau_sd = NULL, tau_rho = NULL, verbose = TRUE ) 
        #- add to final list:
        mcmc_args <- c( mcmc_args, 
           default_list[ !( names( default_list ) %in% names( mcmc_args ) ) ] )
  #}

	res <- list( use_rcpp = use_rcpp, # use rcpp?
		verbose = verbose,
		with_reml = with_reml,
		random_group = random_group,
		fixed_group = fixed_group,
		large = large,
		type_ses = type_ses,
		method = method,
		nlminb_ctrl = nlminb_ctrl,
		with_ses = with_ses,
		raneffs_type = raneffs_type,
		no_pvs = no_pvs,
		mcmc_args = mcmc_args )
	return( res )
}