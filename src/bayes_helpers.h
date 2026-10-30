
#pragma once

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "dmvnorm.h"

Eigen::MatrixXd tsrm_bayes_compute_Vepsi( 
    const Eigen::MatrixXd& Zti,
    const Eigen::MatrixXd& SIGMA_T,
    int nt );

Rcpp::List tsrm_bayes_update_beta_rcpp(
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Zp,
    const Eigen::MatrixXd& Zd,  
    const Eigen::MatrixXd& Zt,
    int nt,
    const Eigen::MatrixXd& G_T,  
    const Eigen::MatrixXd& L_T,
    const Eigen::MatrixXi& groupinfo,
    int ngroups,
    const Rcpp::List& person,  // Liste von Vektoren
    const Rcpp::List& dyads,   // Liste von Vektoren
    const Eigen::MatrixXd& mb,
    const Eigen::MatrixXd& invMb,
    bool large ); 

Rcpp::List tsrm_bayes_update_effects_rcpp(
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Zp,
    const Eigen::MatrixXd& Zd,
    const Eigen::MatrixXd& invG_P,
    const Eigen::MatrixXd& invG_D,
    const Eigen::VectorXd& BETA,
    const Eigen::MatrixXi& groupinfo,
    int ngroups,
    const Rcpp::List& person,  // Liste von Vektoren
    const Rcpp::List& dyads,   // Liste von Vektoren
    SEXP OBJ_invVepsi );