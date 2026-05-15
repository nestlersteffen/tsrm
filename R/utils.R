

srm_include_free_parameters <- function( parm = NULL, 
  parm_list = NULL, parm_table = NULL )
{
  NP <- max(parm_table$index)
  NOP <- nrow(parm_table) 
  for (nn in 1:NOP) {
    free_nn  <- parm_table[nn,]
    type     <- free_nn$type
    pos      <- c( free_nn$pos1, free_nn$pos2 )
    index_nn <- free_nn$index
    x_nn     <- parm_list[[ type ]]   
    x_nn[ pos[1] , pos[2] ] <- parm[ index_nn ]
    if ( type %in% "SD_D" ) {
      x_nn[ pos[1] + 1, pos[2] + 1 ] <- parm[ index_nn ]
    }
    if ( type %in% c("RHO_G", "RHO_P" ) ) {
      x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ]
    }
    if ( type %in% c("RHO_D" ) ) {
      x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ]
      if ( pos[1] == 1 & pos[2] == 3 ) {
        x_nn[ pos[1] + 1 , pos[2] + 1 ] <- parm[ index_nn ]
        x_nn[ pos[2] + 1 , pos[1] + 1 ] <- parm[ index_nn ]
      }
      if ( pos[1] == 1 & pos[2] == 4 ) {
        x_nn[ pos[1] + 1 , pos[2] - 1 ] <- parm[ index_nn ]
        x_nn[ pos[2] - 1 , pos[1] + 1 ] <- parm[ index_nn ]
      }
    }    
    parm_list[[ type ]] <- x_nn       
  }
  return(parm_list)
}

srm_dmvnorm <- function( y = NULL, invSigma = NULL, use_log = TRUE )
{
  p <- dim( invSigma )[1]
  logdetSigma <- -determinant( invSigma, logarithm = TRUE )$modulus[1]
  quadval <- colSums( y * ( invSigma %*% y ) )
  l1 <- - p * log(2*pi) - quadval - logdetSigma
  ll <- 0.5 * l1
  # ... 
  if ( !use_log ) { 
    ll <- exp(ll) 
  }
  return( ll )
}

tsrm_MakeOptFun <- function( grad_fn, ... ) 
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