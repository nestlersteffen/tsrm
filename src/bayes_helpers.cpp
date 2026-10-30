
#include "bayes_helpers.h"

// Eigen::MatrixXd tsrm_bayes_compute_Vepsi( 
//     const Eigen::MatrixXd& Zti,
//     const Eigen::MatrixXd& SIGMA_T,
//     int nt )
// {

//     //- get some information:
//     int p = SIGMA_T.rows();
//     int ni = Zti.rows();
//     //- define output:
//     Eigen::MatrixXd Vepsi = Eigen::MatrixXd::Zero( ni, ni );
//     for ( int k=1; k <= nt; k++ ) {
//         // get submatrix:
//         Eigen::MatrixXd Zk = Zti.block( 0, (k - 1) * p, ni, p );
//         Vepsi += Zk * SIGMA_T * Zk.transpose();
//     }
//     return Vepsi;
// }

Eigen::MatrixXd tsrm_bayes_compute_Vepsi(
    const Eigen::MatrixXd& Zti,
    const Eigen::MatrixXd& L_T, 
    int nt )
{
    const int p  = L_T.rows();
    const int ni = Zti.rows();

    // B = [Z_1 L, Z_2 L, ..., Z_nt L]  (ni x nt*p)
    Eigen::MatrixXd B( ni, nt * p );
    for ( int k = 0; k < nt; ++k ) {
        B.middleCols( k * p, p ).noalias() = Zti.middleCols( k * p, p ) * L_T;
    }

    // Vepsi = B B'  (nur untere Hälfte berechnen, dann spiegeln)
    Eigen::MatrixXd Vepsi = Eigen::MatrixXd::Zero( ni, ni );
    Vepsi.selfadjointView<Eigen::Lower>().rankUpdate( B );
    Vepsi.triangularView<Eigen::StrictlyUpper>() = Vepsi.transpose();
    return Vepsi;
}

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
    bool large ) 
{

    //  std::vector (bleibt im Speicher, für pointer-Positionen)
    std::vector<Eigen::MatrixXd>* invVepsi_vec = new std::vector<Eigen::MatrixXd>(ngroups);

    //- we convert the lists to std::vectors
    std::vector<Eigen::VectorXd> person_vec(ngroups);
    std::vector<Eigen::VectorXd> dyads_vec(ngroups);
    
    for(int ng = 0; ng < ngroups; ng++) {
        person_vec[ng] = Rcpp::as<Eigen::VectorXd>(person[ng]);
        dyads_vec[ng] = Rcpp::as<Eigen::VectorXd>(dyads[ng]);
    }

    //- initialize mu and om:
    Eigen::VectorXd sum_mu = Eigen::VectorXd::Zero(X.cols());
    Eigen::MatrixXd sum_om = Eigen::MatrixXd::Zero(X.cols(), X.cols());
    
    for(int ng = 0; ng < ngroups; ng++) {
            
        // get indices:
        int idx1 = groupinfo(ng, 4) - groupinfo(ng, 3);
        int idx2 = groupinfo(ng, 3);
        
        // get blocks:
        Eigen::VectorXd yi = y.block(idx1, 0, idx2, 1);
        Eigen::MatrixXd Xi = X.block(idx1, 0, idx2, X.cols());
        Eigen::MatrixXd Zpi = Zp.block(idx1, 0, idx2, Zp.cols());
        Eigen::MatrixXd Zdi = Zd.block(idx1, 0, idx2, Zd.cols());
        Eigen::MatrixXd Zti = Zt.block(idx1, 0, idx2, Zt.cols());
            
        // get elements in lists:
        const Eigen::VectorXd& peffi = person_vec[ng];
        const Eigen::VectorXd& deffi = dyads_vec[ng];

        // compute Vepsi:
        Eigen::MatrixXd Vepsi;
        if ( large ) {
            Vepsi = tsrm_bayes_compute_Vepsi( Zti, L_T, nt );
        } else {
            Vepsi = Zti * G_T * Zti.transpose();
        }
        
        // similar to R, we do not use the inverse but the cholesky

        // Eigen::LDLT<Eigen::MatrixXd> ldlt(Vepsi);
        // Eigen::MatrixXd invVepsi = ldlt.solve(Eigen::MatrixXd::Identity(Vepsi.rows(), Vepsi.cols()));
        // (*invVepsi_vec)[ng] = invVepsi;
            
        // Eigen::MatrixXd tmpXi = Xi.transpose() * invVepsi;
        // Eigen::VectorXd residual = yi - Zpi * peffi - Zdi * deffi;
            
        // sum_mu += tmpXi * residual;
        // sum_om += tmpXi * Xi;

        // new version with Eigens backsolve equivalent
        Eigen::LLT<Eigen::MatrixXd> llt(Vepsi);
        Eigen::MatrixXd L = llt.matrixL();          // untere Dreiecksmatrix, V = L L'
        (*invVepsi_vec)[ng] = L;

        // make backsolve variables:
        const int px = Xi.cols();
        Eigen::MatrixXd XR(Xi.rows(), px + 1);
        XR.leftCols(px) = Xi;
        XR.col(px)      = yi - Zpi * peffi - Zdi * deffi;

        // this is backsolve(R, ., transpose = TRUE) )
        L.triangularView<Eigen::Lower>().solveInPlace(XR);

        Eigen::MatrixXd WX = XR.leftCols(px);
        Eigen::VectorXd Wr = XR.col(px);

        sum_mu += WX.transpose() * Wr;              // Xi' V^-1 resid
        sum_om.noalias() += WX.transpose() * WX;               

    }
    
    // compute everything for the BETA update:
    Eigen::MatrixXd OMEGAb = (sum_om + invMb).ldlt().solve( Eigen::MatrixXd::Identity( sum_om.rows(), sum_om.cols() ) );
    Eigen::VectorXd MUb = OMEGAb * sum_mu + OMEGAb * invMb * mb;
    Eigen::VectorXd BETA = rmvnorm_rcpp( MUb, OMEGAb );
    
    // XPtr erstellen (nur 8 Bytes, kein Kopieren!):
    Rcpp::XPtr<std::vector<Eigen::MatrixXd>> ptr(invVepsi_vec, true);
    
    return Rcpp::List::create(
        Rcpp::Named("BETA") = BETA,
        Rcpp::Named("OBJ_invVepsi") = ptr
    );
}

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
    SEXP OBJ_invVepsi )
{

    //- get the saved invVepsi-matrices:
    Rcpp::XPtr<std::vector<Eigen::MatrixXd>> ptr( OBJ_invVepsi );
    std::vector<Eigen::MatrixXd>& invVepsi_vec = *ptr;

    //- we convert the lists to std::vectors
    std::vector<Eigen::VectorXd> person_vec(ngroups);
    std::vector<Eigen::VectorXd> dyads_vec(ngroups);
    
    for(int ng = 0; ng < ngroups; ng++) {
        person_vec[ng] = Rcpp::as<Eigen::VectorXd>(person[ng]);
        dyads_vec[ng] = Rcpp::as<Eigen::VectorXd>(dyads[ng]);
    }

    //- containers to save results
    std::vector<Eigen::VectorXd> new_person_vec(ngroups);
    std::vector<Eigen::VectorXd> new_dyads_vec(ngroups);
    
    // loop through groups:
    for(int ng = 0; ng < ngroups; ng++) {
        
        // get indices:
        int idx1 = groupinfo(ng, 4) - groupinfo(ng, 3);
        int idx2 = groupinfo(ng, 3);
        
        // get blocks:
        Eigen::VectorXd yi = y.block(idx1, 0, idx2, 1);
        Eigen::MatrixXd Xi = X.block(idx1, 0, idx2, X.cols());
        Eigen::MatrixXd Zpi = Zp.block(idx1, 0, idx2, Zp.cols());
        Eigen::MatrixXd Zdi = Zd.block(idx1, 0, idx2, Zd.cols());

        // get elements in lists:
        const Eigen::VectorXd& peffi = person_vec[ng];
        const Eigen::VectorXd& deffi = dyads_vec[ng];

        // make pointer to inverses:
        // const Eigen::MatrixXd& invVepsi = invVepsi_vec[ng];
        
        // Eigen::MatrixXd tmpZVp = Zpi.transpose() * invVepsi;
        // Eigen::MatrixXd tmpZVp_Zpi = tmpZVp * Zpi;
        // Eigen::MatrixXd OMEGAp = ( tmpZVp_Zpi + invG_P ).ldlt().solve(Eigen::MatrixXd::Identity( invG_P.cols(), invG_P.cols()));
        // Eigen::VectorXd residp = yi - Xi*BETA - Zdi * deffi; 
        // Eigen::VectorXd MUp = OMEGAp * ( tmpZVp * residp );
        // Eigen::VectorXd new_peffi = rmvnorm_rcpp( MUp, OMEGAp );

        // Eigen::MatrixXd tmpZVd = Zdi.transpose() * invVepsi;
        // Eigen::MatrixXd tmpZVd_Zdi = tmpZVd * Zdi;
        // Eigen::MatrixXd OMEGAd = ( tmpZVd_Zdi + invG_D ).ldlt().solve(Eigen::MatrixXd::Identity(invG_D.cols(), invG_D.cols()) );
        // Eigen::VectorXd residd = yi - Xi*BETA - Zpi * new_peffi; 
        // Eigen::VectorXd MUd = OMEGAd * (tmpZVd * residd );
        // Eigen::VectorXd new_deffi = rmvnorm_rcpp( MUd, OMEGAd );

        const Eigen::MatrixXd& L = invVepsi_vec[ng];
        const int pzp = Zpi.cols();
        const int pzd = Zdi.cols();
        Eigen::MatrixXd W( Zpi.rows(), pzp + pzd + 1 );
        
        W.leftCols(pzp) = Zpi;
        W.middleCols(pzp,pzd) = Zdi;
        W.col(pzp+pzd) = yi - Xi * BETA;

        L.triangularView<Eigen::Lower>().solveInPlace(W);
        Eigen::MatrixXd Wp = W.leftCols(pzp);
        Eigen::MatrixXd Wd = W.middleCols(pzp,pzd);
        Eigen::VectorXd Wy = W.col(pzp+pzd);

        Eigen::MatrixXd OMEGAp = ( Wp.transpose()*Wp + invG_P ).ldlt().solve(Eigen::MatrixXd::Identity( invG_P.cols(), invG_P.cols()));
        Eigen::VectorXd MUp = OMEGAp * ( Wp.transpose() * ( Wy - Wd*deffi ) );
        Eigen::VectorXd new_peffi = rmvnorm_rcpp( MUp, OMEGAp );

        Eigen::MatrixXd OMEGAd = ( Wd.transpose()*Wd + invG_D ).ldlt().solve(Eigen::MatrixXd::Identity( invG_D.cols(), invG_D.cols()));
        Eigen::VectorXd MUd = OMEGAd * ( Wd.transpose() * ( Wy - Wp*new_peffi ) );
        Eigen::VectorXd new_deffi = rmvnorm_rcpp( MUd, OMEGAd );

        // save result:
        new_person_vec[ng] = new_peffi;
        new_dyads_vec[ng] = new_deffi;

    }

    //- we transform result back to lists
    Rcpp::List new_person(ngroups);
    Rcpp::List new_dyads(ngroups);
    
    for(int ng = 0; ng < ngroups; ng++) {
        new_person[ng] = new_person_vec[ng];
        new_dyads[ng] = new_dyads_vec[ng];
    }
    
    return Rcpp::List::create(
        Rcpp::Named("person") = new_person,
        Rcpp::Named("dyads") = new_dyads
    );
}