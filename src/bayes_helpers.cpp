
#include "bayes_helpers.h"

Rcpp::List tsrm_bayes_update_beta_rcpp(
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Zp,
    const Eigen::MatrixXd& Zd,  
    const Eigen::MatrixXd& Zt,
    const Eigen::MatrixXd& G_T,  
    const Eigen::MatrixXi& groupinfo,
    int ngroups,
    const Rcpp::List& person,  // Liste von Vektoren
    const Rcpp::List& dyads,   // Liste von Vektoren
    const Eigen::MatrixXd& mb,
    const Eigen::MatrixXd& invMb ) 
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
        Eigen::MatrixXd Vepsi = Zti * G_T * Zti.transpose();
            
        // compute inverse and save it in athe pointer vector:
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Vepsi);
        Eigen::MatrixXd invVepsi = ldlt.solve(Eigen::MatrixXd::Identity(Vepsi.rows(), Vepsi.cols()));
        (*invVepsi_vec)[ng] = invVepsi;
            
        // now the sums:
        Eigen::MatrixXd tmpXi = Xi.transpose() * invVepsi;
        Eigen::VectorXd residual = yi - Zpi * peffi - Zdi * deffi;
            
        sum_mu += tmpXi * residual;
        sum_om += tmpXi * Xi;
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
        const Eigen::MatrixXd& invVepsi = invVepsi_vec[ng];
        
        // update person effects:
        Eigen::MatrixXd tmpZVp = Zpi.transpose() * invVepsi;
        Eigen::MatrixXd tmpZVp_Zpi = tmpZVp * Zpi;
        Eigen::MatrixXd OMEGAp = ( tmpZVp_Zpi + invG_P ).ldlt().solve(Eigen::MatrixXd::Identity( invG_P.cols(), invG_P.cols()));
        Eigen::VectorXd residp = yi - Xi*BETA - Zdi * deffi; 
        Eigen::VectorXd MUp = OMEGAp * ( tmpZVp * residp );
        Eigen::VectorXd new_peffi = rmvnorm_rcpp( MUp, OMEGAp );

        // update dyad effects:
        Eigen::MatrixXd tmpZVd = Zdi.transpose() * invVepsi;
        Eigen::MatrixXd tmpZVd_Zdi = tmpZVd * Zdi;
        Eigen::MatrixXd OMEGAd = ( tmpZVd_Zdi + invG_D ).ldlt().solve(Eigen::MatrixXd::Identity(invG_D.cols(), invG_D.cols()) );
        Eigen::VectorXd residd = yi - Xi*BETA - Zpi * new_peffi; 
        Eigen::VectorXd MUd = OMEGAd * (tmpZVd * residd );
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