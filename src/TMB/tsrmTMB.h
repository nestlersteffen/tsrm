namespace tsrmTMB{

    // Matrix Qd for dyadic data: 
    template<class Type>
    matrix<Type> getQd() {
        
        // an empty matrix:
        matrix<Type> Q(6, 6);
        Q.setZero(); 

        // fill the matrix:
        Type c2 = Type(1.0) / sqrt(Type(2.0));
        Q(0, 0) = c2;  Q(0, 1) = c2;
        Q(1, 2) = c2;  Q(1, 3) = c2;
        Q(2, 4) = c2;  Q(2, 5) = c2;
        Q(3, 0) = c2;  Q(3, 1) = -c2;
        Q(4, 2) = c2;  Q(4, 3) = -c2;
        Q(5, 4) = c2;  Q(5, 5) = -c2;
        
        return Q;
    }

    // Matrix Qt for triadic data: 
    template<class Type>
    matrix<Type> getQt() {

        // an empty matrix:
        matrix<Type> Q(6, 6);
        Q.setZero(); 
        
        // fill the matrix:
        Type s6 = std::sqrt(Type(6.0));
        Type s3 = std::sqrt(Type(3.0));

        // row 1
        Q(0, 0) = Type(1.0) / s6;  Q(0, 1) = Type(1.0) / s6;  Q(0, 2) = Type(1.0) / s6;
        Q(0, 3) = Type(1.0) / s6;  Q(0, 4) = Type(1.0) / s6;  Q(0, 5) = Type(1.0) / s6;

        // row 2 
        Q(1, 0) = Type(1.0) / s6;  Q(1, 1) = -Type(1.0) / s6; Q(1, 2) = -Type(1.0) / s6;
        Q(1, 3) = Type(1.0) / s6;  Q(1, 4) = Type(1.0) / s6;  Q(1, 5) = -Type(1.0) / s6;

        // row 3
        Q(2, 0) = Type(1.0) / s3;  Q(2, 1) = -Type(1.0) / s3; Q(2, 2) = Type(1.0) / (Type(2.0) * s3);
        Q(2, 3) = -Type(1.0) / (Type(2.0) * s3); Q(2, 4) = -Type(1.0) / (Type(2.0) * s3); Q(2, 5) = Type(1.0) / (Type(2.0) * s3);

        // row 4
        Q(3, 2) = Type(0.5); Q(3, 3) = -Type(0.5); Q(3, 4) = Type(0.5); Q(3, 5) = -Type(0.5);

        // row 5
        Q(4, 2) = Type(0.5); Q(4, 3) = Type(0.5); Q(4, 4) = -Type(0.5); Q(4, 5) = -Type(0.5);

        // row 6:
        Q(5, 0) = Type(1.0) / s3;  Q(5, 1) = Type(1.0) / s3;  Q(5, 2) = -Type(1.0) / (Type(2.0) * s3);
        Q(5, 3) = -Type(1.0) / (Type(2.0) * s3); Q(5, 4) = -Type(1.0) / (Type(2.0) * s3); Q(5, 5) = -Type(1.0) / (Type(2.0) * s3);

        return Q;
    }

    template<class Type>
    matrix<Type> getS( const vector<Type> &sd ) {
    
        int p = sd.size();
        matrix<Type> S( p, p);
        S.setZero();
        S.diagonal() = sd.exp();
        return S;
    }

    template<class Type>
    matrix<Type> getR2( const vector<Type> &rho ) {
    
        int p = 2;
        matrix<Type> R(p,p);
        R.setIdentity();

        Type r12 = tanh( rho(0) );
        R(0,1) = r12; R(1,0) = r12;
        
        return R;
    }

    template<class Type>
    matrix<Type> getR3( const vector<Type> &rho, int version ) {
    
        int p = 3;
        matrix<Type> R(p,p);
        
        if ( version == 0 ) { // correlation matrix ala TMB
            matrix<Type> L(p,p);
            L.setIdentity();
            int i,j,k=0;
            for(i=0;i<p;i++){
                for(j=0;j<p;j++){
                    if(i>j){L(i,j)=rho(k);k++;}
                }
            }
            matrix<Type> llt=L*L.transpose();
            matrix<Type> tmp=llt;
            for(i=0;i<p;i++){
                for(j=0;j<p;j++){
                    R(i,j)=tmp(i,j)/sqrt(llt(i,i)*llt(j,j));
                }
            }    
        } else if ( version == 1 ) { // a 3 x 3 correlation matrix 
            Type r12 = tanh( rho(0) );
            Type r13 = tanh( rho(1) );
            Type p23 = tanh( rho(2) ); 
            Type r23 = r12*r13 + p23*sqrt( (1-r12*r12)*(1-r13*r13) );

            R.setIdentity();
            R(0,1) = r12; R(1,0) = r12;
            R(0,2) = r13; R(2,0) = r13;
            R(1,2) = r23; R(2,1) = r23;
        }  

        return R;
    }

    // covariance matrix for triadic effects - univariate case:
    template<class Type>
    matrix<Type> uTAU_T( const vector<Type> &sd, const vector<Type> &rho  ) {

        // define matrices
        matrix<Type> S(6,6);
        S.setZero();
        matrix<Type> R(6,6);
        R.setIdentity();
    
        // fill matrix with standard deviations:
        for(int i = 0; i < 6; i++){
          S(i, i) = sd(0);
        }

        // fill matrix with correlations:
        R(0,1) = rho(0); R(1,0) = rho(0);
        R(2,3) = rho(0); R(3,2) = rho(0);
        R(4,5) = rho(0); R(5,4) = rho(0);
        R(0,5) = rho(1); R(5,0) = rho(1);
        R(1,3) = rho(1); R(3,1) = rho(1);
        R(2,4) = rho(1); R(4,2) = rho(1);
        R(0,2) = rho(2); R(2,0) = rho(2);
        R(1,4) = rho(2); R(4,1) = rho(2);
        R(3,5) = rho(2); R(5,3) = rho(2);
        R(0,3) = rho(3); R(3,0) = rho(3);
        R(0,4) = rho(3); R(4,0) = rho(3);
        R(1,2) = rho(3); R(2,1) = rho(3);
        R(1,5) = rho(3); R(5,1) = rho(3);
        R(2,5) = rho(3); R(5,2) = rho(3);
        R(3,4) = rho(3); R(4,3) = rho(3);

        // compute covariance matrix:
        matrix<Type> mat = S * R * S;
        return mat;

    }
        
}