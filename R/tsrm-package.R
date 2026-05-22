#' tsrm: Fit triadic and dyadic round-robin data
#'
#' @description Provides functionality for maximum likelihood, 
#'   restrictced maximum likelihood and Bayesian estimation of 
#'   Social Relations Model Parameters for Triadic and Dyadic Data.
#'
#' @keywords internal
#' @useDynLib tsrm, .registration = TRUE
#' @useDynLib tsrm_TMBExports, .registration = TRUE
#' @importFrom TMB MakeADFun
#' @importFrom methods as
#' @import Rcpp
#' @import RcppEigen
#' @import stats
#' @import utils
#' @import mvtnorm
#' @import numDeriv
#' @import nloptr
"_PACKAGE"