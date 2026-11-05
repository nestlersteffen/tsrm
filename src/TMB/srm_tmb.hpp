/// @file srm_tmb.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type srm_tmb( objective_function<Type>* obj )
{
	
    using namespace density;
	
    // data:
	DATA_MATRIX(y);
    DATA_MATRIX(X);
    DATA_MATRIX(Zp); 
	DATA_MATRIX(Zd);
	DATA_INTEGER(np); 
    DATA_INTEGER(nd);
    DATA_IMATRIX(groupinfo);
    DATA_INTEGER(method); 
  
    // information that we need later:
    int ncolX  = X.cols();
    int ncolZp = Zp.cols();
    int ncolZd = Zd.cols();
	
    // parameters:
    PARAMETER_VECTOR(BETA);
    PARAMETER_VECTOR(vecsdP);
    PARAMETER_VECTOR(vecrhoP);
    PARAMETER_VECTOR(vecsdD);
    PARAMETER_VECTOR(vecrhoD);
  
    // make a covariance matrix for the persons:
    matrix<Type> S_P = tsrmTMB::getS( vecsdP );
    matrix<Type> R_P = tsrmTMB::getR2( vecrhoP );
    matrix<Type> TAU_P = S_P*R_P*S_P;

    // make a covariance matrix for the dyads:
    vector<Type> sdD = vecsdD.exp();
    matrix<Type> S_D = sdD(0) * matrix<Type>::Identity(2,2);
    matrix<Type> R_D = tsrmTMB::getR2( vecrhoD );
    matrix<Type> TAU_D = S_D*R_D*S_D;

    // compute big covariance matrices:
    matrix<Type> identp(np,np);
    identp.setIdentity();
    matrix<Type> SIGMA_P = kronecker(identp,TAU_P);
    
    matrix<Type> identd(nd,nd);
    identd.setIdentity();
    matrix<Type> SIGMA_D = kronecker(identd,TAU_D);
 
    // get group_infos:
    int ngroups = groupinfo.rows();
  
    // iterate through groups to overall likelihood;
    Type nll = 0;
    matrix<Type> reml_penterm(ncolX,ncolX);
    reml_penterm.setZero();
    int idx1 = 0;
    int idx2 = 0;
    for ( int ng = 0; ng < ngroups; ng++ ) {
        
        // how many rows to go:
        idx2 = groupinfo(ng,3);
        
        // now we obtain the smaller matrices:
        vector<Type> ty  = y.block( idx1, 0, idx2, 1 );
        matrix<Type> tX  = X.block( idx1, 0, idx2, ncolX );
        matrix<Type> tZp = Zp.block( idx1, 0, idx2, ncolZp );
        matrix<Type> tZd = Zd.block( idx1, 0, idx2, ncolZd );
        
        // compute observed V:
        matrix<Type> V = tZp*SIGMA_P*tZp.transpose() + tZd*SIGMA_D*tZd.transpose();
        
        // compute mean structure:
        vector<Type> resid = ty - tX*BETA;
    
        // compute normal density:
        nll += MVNORM(V)(resid); 

        // add REML:
        if ( method != 0 ) {
            matrix<Type> iV = V.inverse();
            reml_penterm += tX.transpose()*iV*tX;
        }
    
        // adpapt index:
        idx1 += groupinfo(ng,3);   
    } 
    
    if ( method != 0 ) nll += 0.5*atomic::logdet( reml_penterm );

    // reporting - standard errors:
    vector<Type> outsdP = S_P.diagonal();
    vector<Type> outrhoP(1);
    outrhoP(0) = R_P(0,1); 
    vector<Type> outrhoD(1);
    outrhoD(0) = R_D(0,1);
    
    ADREPORT(BETA);
    ADREPORT(outsdP);
    ADREPORT(outrhoP);
    ADREPORT(sdD);
    ADREPORT(outrhoD);
	
    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this