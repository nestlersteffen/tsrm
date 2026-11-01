
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

Rcpp::List gradient_singlegroup_rcpp(
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
    DeriveFun derivative_fun;
    if ( model == "srm" ) {
        derivative_fun = srm_sigma_derivatives_rcpp;
    } else if ( model == "tsrm" ) {
        derivative_fun = tsrm_sigma_derivatives_rcpp;
    } else if ( model == "htsrm" ) {
        derivative_fun = htsrm_sigma_derivatives_rcpp;
    } else {
        Rcpp::stop("Unknown model '%s' in derivative dispatch.", model);
    }

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
    double llfct = dmvnorm_centered_rcpp( tey, V, true );

    //-- Step 5: derivative for beta:
    Eigen::VectorXd dBETA = tei.transpose() * tX;

    //-- precompute once per group (für jeden Z-Block):
    Eigen::MatrixXd Mg; Eigen::VectorXd wg;
    if ( tZg.cols() > 0 ) {
        Mg = tZg.transpose() * (iV * tZg);
        wg = tZg.transpose() * tei;
    }
    Eigen::MatrixXd Mp = tZp.transpose() * (iV * tZp);
    Eigen::VectorXd wp = tZp.transpose() * tei;
    Eigen::MatrixXd Md = tZd.transpose() * (iV * tZd);
    Eigen::VectorXd wd = tZd.transpose() * tei;
    Eigen::MatrixXd Mt; Eigen::VectorXd wt;
    if ( tZt.cols() > 0 ) {
        Mt = tZt.transpose() * (iV * tZt);
        wt = tZt.transpose() * tei;
    }

    //-- Step 6: REML-Penalty:
    Eigen::MatrixXd penvalue_ll = Eigen::MatrixXd::Zero( tX.cols(), tX.cols() );
    Eigen::MatrixXd XiVZg;
    Eigen::MatrixXd XiVZp;
    Eigen::MatrixXd XiVZd;
    Eigen::MatrixXd XiVZt;
    if ( with_reml ) {
        Eigen::MatrixXd iV_X = ldlt.solve(tX);
        Eigen::MatrixXd iV_Xt = iV_X.transpose();
        Eigen::MatrixXd tXt = tX.transpose();
        penvalue_ll = tXt * iV_X;
        if ( tZg.cols() > 0 ) XiVZg = iV_Xt*tZg;
        XiVZp = iV_Xt*tZp; XiVZd = iV_Xt*tZd;
        if ( tZt.cols() > 0 ) XiVZt = iV_Xt*tZt;
    }
    
    // -- Step 7: iterate across posmat entries for the gradient:
    int NOP = parm_mat.rows();
    int NP  = parm_mat.col(3).maxCoeff();
    Eigen::VectorXd grad = Eigen::VectorXd::Zero( NP ); 
    std::vector<Eigen::MatrixXd> penvalue_gr;
    if ( with_reml ) penvalue_gr.assign( NP, Eigen::MatrixXd::Zero( tX.cols(), tX.cols() ) );
    
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
            Eigen::MatrixXd M;
            Eigen::VectorXd w;
            Eigen::MatrixXd WZ;
                
            if ( type == 1 || type == 2 ) { 
                tmp = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(np, np), sigma_derive );
                tmp = Pp * tmp * Pp.transpose();
                M = Mp; w = wp; WZ = XiVZp;
            } else if ( type == 3 || type == 4 ) {
                tmp = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(nd, nd), sigma_derive );
                tmp = Pd * tmp * Pd.transpose();
                M = Md; w = wd; WZ = XiVZd;
            } else if ( type == 7 || type == 8 ) {
                tmp = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(nt, nt), sigma_derive );
                tmp = Pt * tmp * Pt.transpose();
                M = Mt; w = wt; WZ = XiVZt;
            } else if ( type == 5 || type == 6 ) {
                tmp = sigma_derive;
                M = Mg; w = wg; WZ = XiVZg;
            }

            double pt1 = tmp.cwiseProduct( M ).sum();   // tr(P V') = <S, Z^T P Z>
            double pt2 = w.dot( tmp * w );              // tei^T V' tei = w^T S w
            res = -0.5*( pt1 - pt2 );

            //- things for REML penalty term:
            if ( with_reml ) {
                penvalue_gr[ index - 1] += WZ * tmp * WZ.transpose();
            }

        }

        grad( index - 1 ) += res;

    }
    
    //-- return output:
    return Rcpp::List::create(
        Rcpp::Named("llfct") = llfct,
        Rcpp::Named("gradient") = grad,
        Rcpp::Named("penvalue_ll") = penvalue_ll,
        Rcpp::Named("penvalue_gr") = penvalue_gr
    );
}

Rcpp::List gradient_singlegroup_sparse_rcpp(
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
    DeriveFun derivative_fun;
    if ( model == "srm" ) {
        derivative_fun = srm_sigma_derivatives_rcpp;
    } else if ( model == "tsrm" ) {
        derivative_fun = tsrm_sigma_derivatives_rcpp;
    } else if ( model == "htsrm" ) {
        derivative_fun = htsrm_sigma_derivatives_rcpp;
    } else {
        Rcpp::stop("Unknown model '%s' in derivative dispatch.", model);
    }

    //-- Step 1: compute V:
    Eigen::MatrixXd V = compute_V_sparse(
        tZg, tZp, tZd, tZt,
        SIGMA_G, SIGMA_P, SIGMA_D, SIGMA_T );

    //-- Step 2: compute inverse of V:
    Eigen::LLT<Eigen::MatrixXd> llt(V);
    Eigen::MatrixXd iV = llt.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));

    //-- Step 3: compute "residual":
    Eigen::VectorXd tey = ty - tX * BETA;
    Eigen::VectorXd tei = llt.solve( tey );
    
    //-- Step 4: we also compute the llfct-value:
    double llfct = dmvnorm_centered_rcpp( tey, V, true );
    
    //-- Step 5: derivative for beta:
    Eigen::VectorXd dBETA = tei.transpose() * tX;

    //-- precompute once per group (für jeden Z-Block):
    Eigen::MatrixXd Mg; Eigen::VectorXd wg;
    if ( tZg.cols() > 0 ) {
        Mg = tZg.transpose() * (iV * tZg);
        wg = tZg.transpose() * tei;
    }
    Eigen::MatrixXd Mp = tZp.transpose() * (iV * tZp);
    Eigen::VectorXd wp = tZp.transpose() * tei;
    Eigen::MatrixXd Md = tZd.transpose() * (iV * tZd);
    Eigen::VectorXd wd = tZd.transpose() * tei;
    Eigen::MatrixXd Mt; Eigen::VectorXd wt;
    if ( tZt.cols() > 0 ) {
        Mt = tZt.transpose() * (iV * tZt);
        wt = tZt.transpose() * tei;
    }

    //-- Step 6: REML-Penalty:
    Eigen::MatrixXd penvalue_ll = Eigen::MatrixXd::Zero( tX.cols(), tX.cols() );
    Eigen::MatrixXd XiVZg;
    Eigen::MatrixXd XiVZp;
    Eigen::MatrixXd XiVZd;
    Eigen::MatrixXd XiVZt;
    if ( with_reml ) {
        Eigen::MatrixXd iV_X = llt.solve(tX);
        Eigen::MatrixXd iV_Xt = iV_X.transpose();
        Eigen::MatrixXd tXt = tX.transpose();
        penvalue_ll = tXt * iV_X;
        if ( tZg.cols() > 0 ) XiVZg = iV_Xt*tZg;
        XiVZp = iV_Xt*tZp; XiVZd = iV_Xt*tZd;
        if ( tZt.cols() > 0 ) XiVZt = iV_Xt*tZt;
    }

    // -- Step 7: iterate across posmat entries for the gradient:
    int NOP = parm_mat.rows();
    int NP  = parm_mat.col(3).maxCoeff();
    Eigen::VectorXd grad = Eigen::VectorXd::Zero( NP + 1 ); 
    std::vector<Eigen::MatrixXd> penvalue_gr;
    if ( with_reml ) penvalue_gr.assign( NP, Eigen::MatrixXd::Zero( tX.cols(), tX.cols() ) );

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
            Eigen::MatrixXd tmp;
            Eigen::MatrixXd M;
            Eigen::VectorXd w;
            Eigen::MatrixXd WZ;
                
            if ( type == 1 || type == 2 ) { 
                Eigen::MatrixXd tmp_s = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(np, np), sigma_derive );
                tmp = tmp_s(perm_p,perm_p);
                M = Mp; w = wp; WZ = XiVZp;
            } else if ( type == 3 || type == 4 ) {
                Eigen::MatrixXd tmp_s = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(nd, nd), sigma_derive );
                tmp = tmp_s(perm_d,perm_d);
                M = Md; w = wd; WZ = XiVZd;
            } else if ( type == 7 || type == 8 ) {
                Eigen::MatrixXd tmp_s = Eigen::kroneckerProduct( Eigen::MatrixXd::Identity(nt, nt), sigma_derive );
                tmp = tmp_s(perm_t,perm_t);
                M = Mt; w = wt; WZ = XiVZt;
            } else if ( type == 5 || type == 6 ) {
                tmp = sigma_derive;
                M = Mg; w = wg; WZ = XiVZg;
            }

            // compute gradient value:
            double pt1 = tmp.cwiseProduct( M ).sum();   // tr(P V') = <S, Z^T P Z>
            double pt2 = w.dot( tmp * w );
            res = -0.5*( pt1 - pt2 );
            
            //- things for REML penalty term:
            if ( with_reml ) {
                penvalue_gr[ index - 1] += WZ * tmp * WZ.transpose();
            }

        }

        grad( index - 1 ) += res;

    }
    
    //-- return output:
    return Rcpp::List::create(
        Rcpp::Named("llfct") = llfct,
        Rcpp::Named("gradient") = grad,
        Rcpp::Named("penvalue_ll") = penvalue_ll,
        Rcpp::Named("penvalue_gr") = penvalue_gr
    );
}