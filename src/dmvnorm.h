
#pragma once

//// File Name: dmvnorm.h
//// File Version: 0.01

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <thread>
#include <random>

Eigen::VectorXd rmvnorm_rcpp( const Eigen::VectorXd& mean, 
    const Eigen::MatrixXd& cov );

double log_determinant_rcpp( const Eigen::MatrixXd& M );

double dmvnorm_rcpp( const Eigen::VectorXd& x, 
	                 const Eigen::VectorXd& mu, 
	                 const Eigen::MatrixXd& sigma,
	                 bool use_log = true );

double dmvnorm_centered_rcpp( const Eigen::VectorXd& x, 
                              const Eigen::MatrixXd& sigma,
                              bool use_log = true );