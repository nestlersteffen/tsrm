
#pragma once

//// File Name: derivatives.h
//// File Version: 0.01

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

Eigen::MatrixXd srm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T, 
    const int& type, const Eigen::VectorXi& pos );

Eigen::MatrixXd tsrm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T, 
    const int& type, const Eigen::VectorXi& pos );

Eigen::MatrixXd htsrm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T, 
    const int& type, const Eigen::VectorXi& pos );