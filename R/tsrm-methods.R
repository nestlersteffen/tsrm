
#' Extract Model Coefficients from a glmmX Object
#'
#' Returns the estimated coefficients of a fitted TSRM or SRM.
#'
#' @param object A fitted model object of class `"glmmX"`.
#' @param ... Currently not used.
#'
#' @return A named numeric vector of all model estimates.
#'
#' @seealso [glmm()], [vcov.tsrm()], [summary.tsrm()]
#' @export

coef.tsrm  <- function(object, ...)
{
  return( object$parm )
}

#' Summary method for tsrm objects
#'
#' @param object A fitted model object of class "tsrm"
#' @param digits Number of digits to display
#' @param ... Further arguments passed to or from other methods
#'
#' @seealso [glmm()], [coef.tsrm()], [vcov.tsrm()]
#' @export

summary.tsrm <- function( object, digits = 3L, ...)
{
  output( object=object, digits=digits )
}

#' Variance-Covariance Matrix for a tsrm Object
#'
#' Returns the variance-covariance matrix of the estimated 
#' coefficients of a fitted TSRM or SRM.
#'
#' @param object A fitted model object of class `"tsrm"`.
#' @param ... Currently not used.
#'
#' @return A symmetric numeric matrix. Row and column names match 
#'   [coef.tsrm()].
#'
#' @seealso [glmm()], [coef.tsrm()], [summary.tsrm()]
#' @export

vcov.tsrm <- function(object, ...)
{
  return( object$vcov )
}