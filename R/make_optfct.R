
#---- this function is a wrapper for the optimization routine:

make_optfct <- function( grad_fn, ... ) 
{
	#- the cache
	last_par <- NULL
	last_val <- NULL

	eval_cached <- function(par) {
		if ( !identical(par, last_par) ) {
			last_par <<- par
			last_val <<- grad_fn(par, ... )
		}
		last_val
	}

	list(
		fn = function(par) eval_cached(par)$objective,
		gr = function(par) eval_cached(par)$gradient
	)
}