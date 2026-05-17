
//// File Name: tsrm_exports.cpp
//// File Version: 0.01

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "ml_helpers.h"

// [[Rcpp::export]]
Eigen::MatrixXd srm_compute_V_export(
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    bool random_group )
{
    return srm_compute_V_rcpp( tZg, tZp, tZd, SIGMA_G, SIGMA_P, SIGMA_D );
}

// [[Rcpp::export]]
Eigen::MatrixXd srm_compute_V_sparse_export(
    const Eigen::MatrixXd& tZg,              
    const Eigen::SparseMatrix<double>& tZp,  
    const Eigen::SparseMatrix<double>& tZd,  
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    bool random_group )
{
    return srm_compute_V_sparse( tZg, tZp, tZd, SIGMA_G, SIGMA_P, SIGMA_D );
}

// [[Rcpp::export]]
Eigen::VectorXd srm_gradient_singlegroup_export(
    const Rcpp::List& parm_list,
    const Eigen::MatrixXd& ty,
    const Eigen::MatrixXd& tX,
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::VectorXd& BETA,
    int np, int nd,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml = false, bool random_group = false )
{
    //-- get the parm_list matrices:
    Eigen::MatrixXd SD_G  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_G"]);
    Eigen::MatrixXd SD_P  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_P"]);
    Eigen::MatrixXd SD_D  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_D"]);
    Eigen::MatrixXd RHO_G = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_G"]);
    Eigen::MatrixXd RHO_P = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_P"]);
    Eigen::MatrixXd RHO_D = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_D"]);
    
    //- get result:
    return srm_gradient_singlegroup_rcpp( 
        SD_G, SD_P, SD_D, RHO_G, RHO_P, RHO_D, 
        ty, tX, tZg, tZp, tZd, 
        SIGMA_G, SIGMA_P, SIGMA_D, BETA, 
        np, nd, parm_mat, with_reml, random_group );
}

// [[Rcpp::export]]
Eigen::VectorXd srm_gradient_singlegroup_sparse_export(
    const Rcpp::List& parm_list,
    const Eigen::MatrixXd& ty,
    const Eigen::MatrixXd& tX,
    const Eigen::MatrixXd& tZg,              
    const Eigen::SparseMatrix<double>& tZp,  
    const Eigen::SparseMatrix<double>& tZd,  
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::VectorXd& BETA,
    int np, int nd,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml = false, bool random_group = false )
{
    //-- get the parm_list matrices:
    Eigen::MatrixXd SD_G  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_G"]);
    Eigen::MatrixXd SD_P  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_P"]);
    Eigen::MatrixXd SD_D  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_D"]);
    Eigen::MatrixXd RHO_G = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_G"]);
    Eigen::MatrixXd RHO_P = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_P"]);
    Eigen::MatrixXd RHO_D = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_D"]);
    
    //- get result:
    return srm_gradient_singlegroup_sparse_rcpp( 
        SD_G, SD_P, SD_D, RHO_G, RHO_P, RHO_D, 
        ty, tX, tZg, tZp, tZd, 
        SIGMA_G, SIGMA_P, SIGMA_D, BETA, 
        np, nd, parm_mat, with_reml, random_group );
}