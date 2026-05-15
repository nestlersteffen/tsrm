
#include "ml_helpers.h"

Eigen::MatrixXd srm_compute_V_rcpp(
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    bool random_group )
{
    Eigen::MatrixXd V = tZp * SIGMA_P * tZp.transpose() 
                      + tZd * SIGMA_D * tZd.transpose();
    //- add random group:
    if ( random_group ) {
        V = V + tZg * SIGMA_G * tZg.transpose();
    }
    return V;
}

Eigen::MatrixXd srm_compute_V_sparse(
    const Eigen::MatrixXd& tZg,              // dense
    const Eigen::SparseMatrix<double>& tZp,  // sparse
    const Eigen::SparseMatrix<double>& tZd,  // sparse
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    bool random_group )
{
    typedef Eigen::SparseMatrix<double> SpMat;
    
    SpMat tZp_t = tZp.transpose();
    SpMat tZd_t = tZd.transpose();
    SpMat part_p = tZp * SIGMA_P * tZp_t;
    SpMat part_d = tZd * SIGMA_D * tZd_t;
    SpMat V_sparse = part_p + part_d;
    Eigen::MatrixXd V = Eigen::MatrixXd(V_sparse);
    if (random_group) {
        Eigen::MatrixXd tZg_t = tZg.transpose();
        V += tZg * SIGMA_G * tZg_t;
    }
    return V;
}

Eigen::VectorXd srm_gradient_singlegroup_rcpp(
    const Eigen::MatrixXd& SD_G,      
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D,
    const Eigen::MatrixXd& RHO_G,
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
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
    const Eigen::VectorXi& typevec,
    const Eigen::MatrixXi& posmat,
    bool with_reml, bool random_group ) 
{
    
    //-- Step 1: compute V:
    Eigen::MatrixXd V = srm_compute_V_rcpp( tZg, tZp, tZd, SIGMA_G, SIGMA_P, SIGMA_D, random_group );
        
    //-- Step 2: compute inverse of V:
    Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
    Eigen::MatrixXd iV = ldlt.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));

    //-- Step 3: compute "residual":
    Eigen::VectorXd tey = ty - tX * BETA;
    Eigen::VectorXd tei = iV * tey;
    
    //-- Step 4: we also compute the llfct-value:
    double ng_ll = dmvnorm_centered_rcpp( tey, V, true );

    //-- Step 5: consider REML-Penalty:
    Eigen::MatrixXd P = iV;
    if ( with_reml ) {
        Eigen::MatrixXd iV_X = ldlt.solve(tX);
        Eigen::MatrixXd tXt = tX.transpose();
        Eigen::MatrixXd XtiVX = tXt * iV_X;
        Eigen::LDLT<Eigen::MatrixXd> ldlt2(XtiVX);
        Eigen::MatrixXd iV_Xt = iV_X.transpose();
        Eigen::MatrixXd correction = iV_X * ldlt2.solve(iV_Xt);
        P -= correction;
        ng_ll -= 0.5*log_determinant_rcpp( XtiVX );
    }

    //-- Step 6: derivative for beta:
    Eigen::VectorXd dBETA = tei.transpose() * tX;
    
    // -- Step 7: iterate across posmat entries for the gradient:
    int no_parm = posmat.rows();
    Eigen::VectorXd ng_grad = Eigen::VectorXd::Zero( no_parm + 1 ); 

    for(int nn = 0; nn < no_parm; nn++) {

        // get type of parameter:
        int type = typevec(nn);
        Eigen::VectorXi pos = posmat.row(nn);

        // make a placeholder:
        double res = 0;

        // now compute the gradient values:
        if ( type == 0 ) {
                
            res = dBETA( pos(0) - 1 );
            
        } else { // one of the covariance matrices:

            //- get the derivative:
            Eigen::MatrixXd sigma_derive = srm_sigma_derivatives_rcpp( 
                SD_G, SD_P, SD_D, RHO_G, RHO_P, RHO_D, type, pos );
    
            // //- compute derivative of V depending on matrix:
            Eigen::MatrixXd tmp;
            Eigen::MatrixXd Z;
                
            if ( type == 1 || type == 2 ) { 
                Eigen::MatrixXd tmp_p( np, np );
                tmp_p.setIdentity();
                tmp = Eigen::kroneckerProduct( tmp_p, sigma_derive );
                Z = tZp;
            } else if ( type == 3 || type == 4 ) {
                Eigen::MatrixXd tmp_d( nd, nd );
                tmp_d.setIdentity();
                tmp = Eigen::kroneckerProduct( tmp_d, sigma_derive );
                Z = tZd;
            } else if ( type == 5 || type == 6 ) { 
                tmp = sigma_derive;
                Z = tZg;
            }

            Eigen::MatrixXd Z_t = Z.transpose();
            Eigen::MatrixXd ZS = Z * tmp;
            Eigen::MatrixXd V_DERIVE = ZS * Z_t;
            Eigen::VectorXd Vtei = V_DERIVE * tei;

            // compute gradient value:
            double pt1 = P.cwiseProduct( V_DERIVE ).sum();
            // double pt2 = (tei.transpose() * V_DERIVE * tei).value();
            double pt2 = tei.dot( Vtei );
            res = -0.5*( pt1 - pt2 );

        }

        ng_grad( nn + 1 ) = res;

    }

    //-- add ll-value:
    ng_grad( 0 ) = ng_ll;
    
    //-- return output:
    return ng_grad;
}

Eigen::VectorXd srm_gradient_singlegroup_sparse_rcpp(
    const Eigen::MatrixXd& SD_G,      
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D,
    const Eigen::MatrixXd& RHO_G,
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
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
    const Eigen::VectorXi& typevec,
    const Eigen::MatrixXi& posmat,
    bool with_reml, bool random_group ) 
{
    
    //-- Step 1: compute V:
    Eigen::MatrixXd V = srm_compute_V_sparse( tZg, tZp, tZd, SIGMA_G, SIGMA_P, SIGMA_D, random_group );

    //-- Step 2: compute inverse of V:
    Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
    Eigen::MatrixXd iV = ldlt.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));

    //-- Step 3: compute "residual":
    Eigen::VectorXd tey = ty - tX * BETA;
    Eigen::VectorXd tei = iV * tey;
    
    //-- Step 4: we also compute the llfct-value:
    double ng_ll = dmvnorm_centered_rcpp( tey, V, true );
    
    //-- Step 5: consider REML-Penalty:
    Eigen::MatrixXd P = iV;
    if ( with_reml ) {
        Eigen::MatrixXd iV_X = ldlt.solve(tX);
        Eigen::MatrixXd tXt = tX.transpose();
        Eigen::MatrixXd XtiVX = tXt * iV_X;
        Eigen::LDLT<Eigen::MatrixXd> ldlt2(XtiVX);
        Eigen::MatrixXd iV_Xt = iV_X.transpose();
        Eigen::MatrixXd correction = iV_X * ldlt2.solve(iV_Xt);
        P -= correction;
        ng_ll -= 0.5*log_determinant_rcpp( XtiVX );
    }

    //-- Step 6: derivative for beta:
    Eigen::VectorXd dBETA = tei.transpose() * tX;
    
    // -- Step 7: iterate across posmat entries for the gradient:
    int no_parm = posmat.rows();
    Eigen::VectorXd ng_grad = Eigen::VectorXd::Zero( no_parm + 1 ); 

    for(int nn = 0; nn < no_parm; nn++) {

        // get type of parameter:
        int type = typevec(nn);
        Eigen::VectorXi pos = posmat.row(nn);

        // make a placeholder:
        double res = 0;

        // now compute the gradient values:
        if ( type == 0 ) {
                
            res = dBETA( pos(0) - 1 );
            
        } else { // one of the covariance matrices:

            //- get the derivative:
            Eigen::MatrixXd sigma_derive = srm_sigma_derivatives_rcpp( 
                SD_G, SD_P, SD_D, RHO_G, RHO_P, RHO_D, type, pos );
    
            //- compute derivative of V depending on matrix:
            Eigen::MatrixXd V_DERIVE;
                
            if ( type == 1 || type == 2 ) { 
                SpMat tmp_p( np, np );
                tmp_p.setIdentity();
                SpMat tmp = Eigen::kroneckerProduct( tmp_p, sigma_derive.sparseView() ).eval();
                SpMat tZp_t = tZp.transpose();
                SpMat ZS = tZp * tmp;
                SpMat V_DERIVE_sparse = ZS * tZp_t;
                V_DERIVE = Eigen::MatrixXd( V_DERIVE_sparse );
            } else if ( type == 3 || type == 4 ) {
                SpMat tmp_d( nd, nd );
                tmp_d.setIdentity();
                SpMat tmp = Eigen::kroneckerProduct( tmp_d, sigma_derive.sparseView() ).eval();
                SpMat tZd_t = tZd.transpose();
                SpMat ZS = tZd * tmp;
                SpMat V_DERIVE_sparse = ZS * tZd_t;
                V_DERIVE = Eigen::MatrixXd( V_DERIVE_sparse );
            } else if ( type == 5 || type == 6 ) {
                Eigen::MatrixXd tZg_t = tZg.transpose();
                Eigen::MatrixXd ZSg = tZg * sigma_derive;
                V_DERIVE = ZSg * tZg_t;
            }

            Eigen::VectorXd Vtei = V_DERIVE * tei;

            // compute gradient value:
            double pt1 = P.cwiseProduct( V_DERIVE ).sum();
            double pt2 = tei.dot( Vtei );
            res = -0.5*( pt1 - pt2 );
            
        }

        ng_grad( nn + 1 ) = res;

    }

    //-- add ll-value:
    ng_grad( 0 ) = ng_ll;
    
    //-- return output:
    return ng_grad;
}