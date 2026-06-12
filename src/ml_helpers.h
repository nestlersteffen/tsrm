
#pragma once

//// File Name: ml_helpers.h
//// File Version: 0.01

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <unsupported/Eigen/KroneckerProduct>
#include "derivatives.h"
#include "dmvnorm.h"

typedef Eigen::SparseMatrix<double> SpMat;

Eigen::MatrixXd compute_V_rcpp(
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& tZt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::MatrixXd& SIGMA_T );

Eigen::MatrixXd compute_V_sparse(
    const Eigen::MatrixXd& tZg,              // dense
    const Eigen::SparseMatrix<double>& tZp,  // sparse
    const Eigen::SparseMatrix<double>& tZd,  // sparse
    const Eigen::SparseMatrix<double>& tZt,  // sparse
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T );

Eigen::VectorXd gradient_singlegroup_rcpp(
    const Eigen::MatrixXd& SD_G,      
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D,
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G,
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T,
    const Eigen::MatrixXd& ty,
    const Eigen::MatrixXd& tX,
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& tZt,
    const Eigen::MatrixXd& Pp, 
    const Eigen::MatrixXd& Pd,
    const Eigen::MatrixXd& Pt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::MatrixXd& SIGMA_T,
    const Eigen::VectorXd& BETA,
    int np, int nd, int nt,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml = false, const std::string& model = "srm" );

Eigen::VectorXd gradient_singlegroup_sparse_rcpp(
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D,
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G,
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T,
    const Eigen::MatrixXd& ty,
    const Eigen::MatrixXd& tX,
    const Eigen::MatrixXd& tZg,
    const Eigen::SparseMatrix<double>& tZp,  
    const Eigen::SparseMatrix<double>& tZd,  
    const Eigen::SparseMatrix<double>& tZt,
    const Eigen::SparseMatrix<double>& Pp, 
    const Eigen::SparseMatrix<double>& Pd,
    const Eigen::SparseMatrix<double>& Pt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T,
    const Eigen::VectorXd& BETA,
    int np, int nd, int nt,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml = false, const std::string& model = "srm" );