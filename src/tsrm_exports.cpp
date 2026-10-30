
//// File Name: tsrm_exports.cpp
//// File Version: 0.01

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "ml_helpers.h"
#include "bayes_helpers.h"

// [[Rcpp::export]]
void tsrm_set_eigen_threads( int n ) {
    Eigen::setNbThreads( n );
}

// [[Rcpp::export]]
int tsrm_get_eigen_threads() {
    return Eigen::nbThreads();
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_V_export(
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& tZt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::MatrixXd& SIGMA_T )
{
    return compute_V_rcpp( tZg, tZp, tZd, tZt, SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T );
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_V_sparse_export(
    const Eigen::MatrixXd& tZg,
    const Eigen::SparseMatrix<double>& tZp,
    const Eigen::SparseMatrix<double>& tZd,
    const Eigen::SparseMatrix<double>& tZt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T )
{
    return compute_V_sparse( tZg, tZp, tZd, tZt, SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T );
}

// [[Rcpp::export]]
Eigen::VectorXd gradient_singlegroup_export(
    const Rcpp::List& parm_list,
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
    bool with_reml, const std::string& model )
{
    
    //-- get SD and RHO matrices from parm_list:
    Eigen::MatrixXd SD_P  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_P"]);
    Eigen::MatrixXd SD_D  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_D"]);
    Eigen::MatrixXd RHO_P = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_P"]);
    Eigen::MatrixXd RHO_D = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_D"]);
    
    //-- model-specific matrices:
    Eigen::MatrixXd SD_G, SD_T, RHO_G, RHO_T;
    if ( model == "srm" ) {
        SD_G  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_G"]);
        RHO_G = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_G"]);
        SD_T  = Eigen::MatrixXd(0, 0);
        RHO_T = Eigen::MatrixXd(0, 0);
    } else {
        SD_T  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_T"]);
        RHO_T = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_T"]);
        SD_G  = Eigen::MatrixXd(0, 0);
        RHO_G = Eigen::MatrixXd(0, 0);
    }

    return gradient_singlegroup_rcpp(
        SD_G, SD_P, SD_D, SD_T,
        RHO_G, RHO_P, RHO_D, RHO_T,
        ty, tX, tZg, tZp, tZd, tZt, 
        Pp, Pd, Pt,
        SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T,
        BETA, np, nd, nt,
        parm_mat, with_reml, model );
}

// [[Rcpp::export]]
Eigen::VectorXd gradient_singlegroup_sparse_export(
    const Rcpp::List& parm_list,
    const Eigen::MatrixXd& ty,
    const Eigen::MatrixXd& tX,
    const Eigen::MatrixXd& tZg,
    const Eigen::SparseMatrix<double>& tZp,
    const Eigen::SparseMatrix<double>& tZd,
    const Eigen::SparseMatrix<double>& tZt,
    const Eigen::MatrixXd& Pp, 
    const Eigen::MatrixXd& Pd,
    const Eigen::MatrixXd& Pt,
    const Eigen::VectorXi& perm_p, 
    const Eigen::VectorXi& perm_d,
    const Eigen::VectorXi& perm_t,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T,
    const Eigen::VectorXd& BETA,
    int np, int nd, int nt,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml,
    const std::string& model )
{
    
    Eigen::MatrixXd SD_P  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_P"]);
    Eigen::MatrixXd SD_D  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_D"]);
    Eigen::MatrixXd RHO_P = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_P"]);
    Eigen::MatrixXd RHO_D = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_D"]);

    Eigen::MatrixXd SD_G, SD_T, RHO_G, RHO_T;
    if ( model == "srm" ) {
        SD_G  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_G"]);
        RHO_G = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_G"]);
        SD_T  = Eigen::MatrixXd(0, 0);
        RHO_T = Eigen::MatrixXd(0, 0);
    } else {
        SD_T  = Rcpp::as<Eigen::MatrixXd>(parm_list["SD_T"]);
        RHO_T = Rcpp::as<Eigen::MatrixXd>(parm_list["RHO_T"]);
        SD_G  = Eigen::MatrixXd(0, 0);
        RHO_G = Eigen::MatrixXd(0, 0);
    }

    return gradient_singlegroup_sparse_rcpp(
        SD_G, SD_P, SD_D, SD_T,
        RHO_G, RHO_P, RHO_D, RHO_T,
        ty, tX, tZg, tZp, tZd, tZt, 
        Pp, Pd, Pt,
        perm_p, perm_d, perm_t,
        SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T,
        BETA, np, nd, nt,
        parm_mat, with_reml, model );
}

// [[Rcpp::export]]
Rcpp::List tsrm_bayes_update_beta_export(
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
    bool large )
{

    return tsrm_bayes_update_beta_rcpp(
        y, X, Zp, Zd, Zt, nt, G_T, L_T, groupinfo,
        ngroups, person, dyads, mb, invMb, large );

} 

// [[Rcpp::export]]
Rcpp::List tsrm_bayes_update_effects_export(
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
    SEXP OBJ_invVepsi )
{

    return tsrm_bayes_update_effects_rcpp( 
        y, X, Zp, Zd, invG_P, invG_D, BETA, 
        groupinfo, ngroups, person, dyads, OBJ_invVepsi );

}