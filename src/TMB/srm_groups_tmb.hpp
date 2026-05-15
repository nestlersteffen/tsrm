/// @file srm_groups_tmb.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type srm_groups_tmb( objective_function<Type>* obj )
{
	
    using namespace density;
  	
    // data:
  	DATA_MATRIX(y);
    DATA_MATRIX(X);
    DATA_MATRIX(Zg);
    DATA_MATRIX(Zp); 
  	DATA_MATRIX(Zd);
  	DATA_INTEGER(np); 
    DATA_INTEGER(nd);
    DATA_IMATRIX(groupinfo);
    DATA_INTEGER(method); 
    DATA_INTEGER(no_var);
    
    // information that we need later:
    int ncolX  = X.cols();
    int ncolZg = Zg.cols();
    int ncolZp = Zp.cols();
    int ncolZd = Zd.cols();
  	
    // parameters:
    PARAMETER_VECTOR(BETA);
    PARAMETER_VECTOR(sdP);
    PARAMETER_VECTOR(rhoP);
    PARAMETER_VECTOR(sdD);
    PARAMETER_VECTOR(rhoD);
    PARAMETER_VECTOR(sdG);
    PARAMETER_VECTOR(rhoG);
    
    // make a covariance matrix for the persons:
    matrix<Type> TAU_P = tsrmTMB::TAU( sdP, rhoP );
    
    // make a covariance matrix for the dyads:
    matrix<Type> TAU_D = tsrmTMB::bTAU_D( sdD, rhoD, no_var );

    // make covariance matrix for the groups:
    matrix<Type> TAU_G = tsrmTMB::TAU( sdG, rhoG );

    // compute big covariance matrices:
    matrix<Type> identp(np,np);
    identp.setIdentity();
    matrix<Type> SIGMA_P = kronecker(TAU_P, identp);
    
    matrix<Type> identd(nd,nd);
    identd.setIdentity();
    matrix<Type> SIGMA_D = kronecker(TAU_D, identd);
   
    // get group_infos:
    int ngroups = groupinfo.rows();
    
    // iterate through groups to overall likelihood;
    Type nll = 0;
    int idx1 = 0;
    int idx2 = 0;
    
    for ( int ng = 0; ng < ngroups; ng++ ) {
        // how many rows to go:
        idx2 = groupinfo(ng,2);
      
        // now we obtain the smaller matrices:
        vector<Type> ty  = y.block( idx1, 0, idx2, 1 );
        matrix<Type> tX  = X.block( idx1, 0, idx2, ncolX );
        matrix<Type> tZg = Zg.block( idx1, 0, idx2, ncolZg );
        matrix<Type> tZp = Zp.block( idx1, 0, idx2, ncolZp );
        matrix<Type> tZd = Zd.block( idx1, 0, idx2, ncolZd );
        
        // compute observed V:
        matrix<Type> V = tZg*TAU_G*tZg.transpose() + tZp*SIGMA_P*tZp.transpose() + tZd*SIGMA_D*tZd.transpose();
        
        // compute mean structure:
        vector<Type> resid = ty - tX*BETA;
        
        // compute normal density:
        nll += MVNORM(V)(resid); 
        
        // add REML:
        if ( method != 0 ) {
            matrix<Type> iV = V.inverse();
            matrix<Type> tmp = tX.transpose()*iV*tX;
            nll += 0.5*atomic::logdet( tmp );
        }
        // adpapt index:
        idx1 += groupinfo(ng,2);   
    } 
    
    // reporting:
    ADREPORT(BETA);
    ADREPORT(TAU_G);
    ADREPORT(TAU_P);
    ADREPORT(TAU_D);
  	
    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this