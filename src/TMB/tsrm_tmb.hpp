/// @file tsrm_tmb.hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tsrm_tmb( objective_function<Type>* obj )
{
	
    using namespace density;
    
    // data:
	DATA_MATRIX(y);
    DATA_MATRIX(X);
    DATA_MATRIX(Zp); 
	DATA_MATRIX(Zd);
    DATA_MATRIX(Zt);
	DATA_INTEGER(np); 
    DATA_INTEGER(nd);
    DATA_INTEGER(nt);
    DATA_IMATRIX(groupinfo);
    DATA_INTEGER(method); 
    DATA_INTEGER(version);  
    
    // information that we need later:
    int ncolX  = X.cols();
    int ncolZp = Zp.cols();
    int ncolZd = Zd.cols();
	int ncolZt = Zt.cols();

    // parameters:
    PARAMETER_VECTOR(BETA);
    PARAMETER_VECTOR(vecsdP);
    PARAMETER_VECTOR(vecrhoP);
    PARAMETER_VECTOR(vecsdD);
    PARAMETER_VECTOR(vecrhoD);
    PARAMETER_VECTOR(vecsdT);
    PARAMETER_VECTOR(vecrhoT);

    // make a covariance matrix for the persons:
    matrix<Type> S_P = tsrmTMB::getS( vecsdP );
    matrix<Type> R_P = tsrmTMB::getR3( vecrhoP, version );
    matrix<Type> TAU_P = S_P*R_P*S_P;
  
    // make a covariance matrix for the dyads:
    vector<Type> vecsdD1 = vecsdD.head(3);
    vector<Type> vecsdD2 = vecsdD.tail(3);
    matrix<Type> S_D1 = tsrmTMB::getS( vecsdD1 );
    matrix<Type> S_D2 = tsrmTMB::getS( vecsdD2 );
    vector<Type> vecrhoD1 = vecrhoD.head(3);
    vector<Type> vecrhoD2 = vecrhoD.tail(3);
    matrix<Type> R_D1 = tsrmTMB::getR3( vecrhoD1, version );
    matrix<Type> R_D2 = tsrmTMB::getR3( vecrhoD2, version );
    matrix<Type> SR_D1 = S_D1*R_D1*S_D1;
    matrix<Type> SR_D2 = S_D2*R_D2*S_D2;

    // make final TAU_D:
    matrix<Type> SR_D(6,6);
    SR_D.setZero(); 
    SR_D.block(0, 0, 3, 3) = SR_D1;        // block(startRow, startCol, nRows, nCols)
    SR_D.block(3, 3, 3, 3) = SR_D2;

    // get Q:
    matrix<Type> Qd = tsrmTMB::getQd<Type>();
    matrix<Type> TAU_D = Qd.transpose()*SR_D*Qd;

    // make a covariance matrix for the triads:
    matrix<Type> TAU_T = tsrmTMB::uTAU_T( vecsdT, vecrhoT );
   
    // compute big covariance matrices:
    matrix<Type> identp(np,np);
    identp.setIdentity();
    matrix<Type> SIGMA_P = kronecker(identp,TAU_P);
    
    matrix<Type> identd(nd,nd);
    identd.setIdentity();
    matrix<Type> SIGMA_D = kronecker(identd,TAU_D);

    matrix<Type> identt(nt,nt);
    identt.setIdentity();
    matrix<Type> SIGMA_T = kronecker(identt,TAU_T);
 
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
        matrix<Type> tZt = Zt.block( idx1, 0, idx2, ncolZt );
        
        // compute observed V:
        matrix<Type> V = tZp*SIGMA_P*tZp.transpose() + tZd*SIGMA_D*tZd.transpose() + tZt*SIGMA_T*tZt.transpose();
        
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
    vector<Type> outrhoP(3);
    outrhoP(0) = R_P(0,1);
    outrhoP(1) = R_P(0,2);
    outrhoP(2) = R_P(1,2); 
    
    vector<Type> sdD_all = TAU_D.diagonal().array().sqrt();
    matrix<Type> out_invS_D(6,6);
    out_invS_D.setZero();
    out_invS_D.diagonal() = sdD_all.inverse();
    vector<Type> outsdD(3);                                   
    outsdD(0) = sdD_all(0);
    outsdD(1) = sdD_all(2);
    outsdD(2) = sdD_all(4);
    matrix<Type> outR_D = out_invS_D*TAU_D*out_invS_D; 
    vector<Type> outrhoD(9);
    outrhoD(0) = outR_D(0,1); 
    outrhoD(1) = outR_D(2,3);
    outrhoD(2) = outR_D(4,5);
    outrhoD(3) = outR_D(0,2);
    outrhoD(4) = outR_D(0,5);
    outrhoD(5) = outR_D(2,4);
    outrhoD(6) = outR_D(0,3);
    outrhoD(7) = outR_D(0,4);
    outrhoD(8) = outR_D(2,5);

    ADREPORT(BETA);
    ADREPORT(outsdP);
    ADREPORT(outrhoP);
    ADREPORT(outsdD);
    ADREPORT(outrhoD);
    ADREPORT(vecsdT);
    ADREPORT(vecrhoT);

    // reporting - checks:
    REPORT(BETA);
    REPORT(outsdP);
    REPORT(outrhoP);
    REPORT(outsdD);
    REPORT(outrhoD);
    REPORT(vecsdT);
    REPORT(vecrhoT);

	return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this