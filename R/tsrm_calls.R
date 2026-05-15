
coef.tsrm  <- function(object, ...)
{
  return( object$parm )
}

summary.tsrm <- function( object, digits = 3L, ...)
{
  srm_out( object = object, digits = digits )
}

vcov.mptmem <- function(object, ...)
{
  return( object$vcov )
}