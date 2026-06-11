
#include "ml_helpers.h"

Eigen::MatrixXd compute_V_rcpp(
    const Eigen::MatrixXd& tZg,
    const Eigen::MatrixXd& tZp,
    const Eigen::MatrixXd& tZd,
    const Eigen::MatrixXd& tZt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::MatrixXd& SIGMA_T )
{
    Eigen::MatrixXd V = tZp * SIGMA_P * tZp.transpose() 
                      + tZd * SIGMA_D * tZd.transpose();
    //- add group or triadic effects:
    if ( tZt.cols() > 0 ) {
        V += tZt * SIGMA_T * tZt.transpose();
    }
    if ( tZg.cols() > 0 ) {
        V += tZg * SIGMA_G * tZg.transpose();
    }
    return V;
}

Eigen::MatrixXd compute_V_sparse(
    const Eigen::MatrixXd& tZg,
    const Eigen::SparseMatrix<double>& tZp,
    const Eigen::SparseMatrix<double>& tZd,
    const Eigen::SparseMatrix<double>& tZt,
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T )
{
    SpMat tZp_t = tZp.transpose();
    SpMat tZd_t = tZd.transpose();
    SpMat part_p = tZp * SIGMA_P * tZp_t;
    SpMat part_d = tZd * SIGMA_D * tZd_t;
    SpMat V_sparse = part_p + part_d;
    
    if ( tZt.cols() > 0 ) {
        SpMat tZt_t = tZt.transpose();
        V_sparse += tZt * SIGMA_T * tZt_t;
    }

    Eigen::MatrixXd V = Eigen::MatrixXd(V_sparse);

    if ( tZg.cols() > 0 ) {
        Eigen::MatrixXd tZg_t = tZg.transpose();
        V += tZg * SIGMA_G * tZg_t;
    }
    return V;
}

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
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::MatrixXd& SIGMA_P,
    const Eigen::MatrixXd& SIGMA_D,
    const Eigen::MatrixXd& SIGMA_T,
    const Eigen::VectorXd& BETA,
    int np, int nd, int nt,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml, const std::string& model ) 
{
    
    //-- set pointer to correct derivatives function:
    typedef Eigen::MatrixXd (*DeriveFun)(
        const Eigen::MatrixXd& SD_G,
        const Eigen::MatrixXd& SD_P,
        const Eigen::MatrixXd& SD_D,
        const Eigen::MatrixXd& SD_T,
        const Eigen::MatrixXd& RHO_G,
        const Eigen::MatrixXd& RHO_P,
        const Eigen::MatrixXd& RHO_D,
        const Eigen::MatrixXd& RHO_T,
        const int& type, const Eigen::VectorXi& pos
    );
    DeriveFun derivative_fun = (model == "srm")
        ? srm_sigma_derivatives_rcpp
        : tsrm_sigma_derivatives_rcpp;

    //-- Step 1: compute V:
    Eigen::MatrixXd V = tZp * SIGMA_P * tZp.transpose()
                      + tZd * SIGMA_D * tZd.transpose();
    if ( tZt.cols() > 0 ) {
        V += tZt * SIGMA_T * tZt.transpose();
    }
    if ( tZg.cols() > 0 ) {
        V += tZg * SIGMA_G * tZg.transpose();
    }
        
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
    int NOP = parm_mat.rows();
    int NP  = parm_mat.col(3).maxCoeff();
    Eigen::VectorXd ng_grad = Eigen::VectorXd::Zero( NP + 1 ); 
    for(int nn = 0; nn < NOP; nn++) {

        // get type of parameter:
        int type  = parm_mat(nn,2);
        int index = parm_mat(nn,3);
        Eigen::VectorXi pos = parm_mat.row(nn).segment<2>(0);

        // make a placeholder:
        double res = 0;

        // now compute the gradient values:
        if ( type == 0 ) {
                
            res = dBETA( pos(0) - 1 );
            
        } else { // one of the covariance matrices:

            //- get the derivative:
            Eigen::MatrixXd sigma_derive = derivative_fun( 
                SD_G, SD_P, SD_D, SD_T,
                RHO_G, RHO_P, RHO_D, RHO_T,
                type, pos );
    
            // //- compute derivative of V depending on matrix:
            Eigen::MatrixXd tmp;
            Eigen::MatrixXd Z;
                
            if ( type == 1 || type == 2 ) {
                tmp = Eigen::kroneckerProduct(
                    Eigen::MatrixXd::Identity(np, np), sigma_derive );
                tmp = Pp * tmp * Pp.transpose();
                Z = tZp;
            } else if ( type == 3 || type == 4 ) {
                tmp = Eigen::kroneckerProduct(
                    Eigen::MatrixXd::Identity(nd, nd), sigma_derive );
                tmp = Pd * tmp * Pd.transpose();
                Z = tZd;
            } else if ( type == 7 || type == 8 ) {
                tmp = Eigen::kroneckerProduct(
                    Eigen::MatrixXd::Identity(nt, nt), sigma_derive );
                Z = tZt;
            } else if ( type == 5 || type == 6 ) {
                tmp = sigma_derive;
                Z   = tZg;
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

        ng_grad( index ) = res;

    }

    //-- add ll-value:
    ng_grad( 0 ) = ng_ll;
    
    //-- return output:
    return ng_grad;
}

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
    const Eigen::MatrixXd& SIGMA_G,
    const Eigen::SparseMatrix<double>& SIGMA_P,
    const Eigen::SparseMatrix<double>& SIGMA_D,
    const Eigen::SparseMatrix<double>& SIGMA_T,
    const Eigen::VectorXd& BETA,
    int np, int nd, int nt,
    const Eigen::MatrixXi& parm_mat,
    bool with_reml, const std::string& model ) 
{
    
    //-- set pointer to correct derivatives function:
    typedef Eigen::MatrixXd (*DeriveFun)(
        const Eigen::MatrixXd& SD_G,
        const Eigen::MatrixXd& SD_P,
        const Eigen::MatrixXd& SD_D,
        const Eigen::MatrixXd& SD_T,
        const Eigen::MatrixXd& RHO_G,
        const Eigen::MatrixXd& RHO_P,
        const Eigen::MatrixXd& RHO_D,
        const Eigen::MatrixXd& RHO_T,
        const int& type, const Eigen::VectorXi& pos
    );
    DeriveFun derivative_fun = (model == "srm")
        ? srm_sigma_derivatives_rcpp
        : tsrm_sigma_derivatives_rcpp;

    //-- Step 1: compute V:
    Eigen::MatrixXd V = compute_V_sparse(
        tZg, tZp, tZd, tZt,
        SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T );

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
    int NOP = parm_mat.rows();
    int NP  = parm_mat.col(3).maxCoeff();
    Eigen::VectorXd ng_grad = Eigen::VectorXd::Zero( NP + 1 ); 
    for(int nn = 0; nn < NOP; nn++) {

        // get type of parameter:
        int type = parm_mat(nn,2);
        int index = parm_mat(nn,3);
        Eigen::VectorXi pos = parm_mat.row(nn).segment<2>(0);

        // make a placeholder:
        double res = 0;

        // now compute the gradient values:
        if ( type == 0 ) {
                
            res = dBETA( pos(0) - 1 );
            
        } else { // one of the covariance matrices:

            //- get the derivative:
            Eigen::MatrixXd sigma_derive = derivative_fun(
                SD_G, SD_P, SD_D, SD_T,
                RHO_G, RHO_P, RHO_D, RHO_T,
                type, pos );
    
            //- compute derivative of V depending on matrix:
            Eigen::MatrixXd V_DERIVE;
                
            if ( type == 1 || type == 2 ) { 
                SpMat tmp_p( np, np );
                tmp_p.setIdentity();
                SpMat tmp = Eigen::kroneckerProduct( tmp_p, sigma_derive.sparseView() ).eval();
                tmp = Pp * tmp * Pp.transpose();
                SpMat tZp_t = tZp.transpose();
                SpMat ZS = tZp * tmp;
                SpMat V_DERIVE_sparse = ZS * tZp_t;
                V_DERIVE = Eigen::MatrixXd( V_DERIVE_sparse );
            } else if ( type == 3 || type == 4 ) {
                SpMat tmp_d( nd, nd );
                tmp_d.setIdentity();
                SpMat tmp = Eigen::kroneckerProduct( tmp_d, sigma_derive.sparseView() ).eval();
                tmp = Pd * tmp * Pd.transpose() ;
                SpMat tZd_t = tZd.transpose();
                SpMat ZS = tZd * tmp;
                SpMat V_DERIVE_sparse = ZS * tZd_t;
                V_DERIVE = Eigen::MatrixXd( V_DERIVE_sparse );
            } else if ( type == 7 || type == 8 ) {
                SpMat tmp_t( nt, nt );
                tmp_t.setIdentity();
                SpMat tmp = Eigen::kroneckerProduct( tmp_t, sigma_derive.sparseView() ).eval();
                SpMat tZt_t = tZt.transpose();
                SpMat ZS = tZt * tmp;
                SpMat V_DERIVE_sparse = ZS * tZt_t;
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

        ng_grad( index ) = res;

    }

    //-- add ll-value:
    ng_grad( 0 ) = ng_ll;
    
    //-- return output:
    return ng_grad;
}