
namespace tsrmTMB{

    // standard covariance matrix given standard deviations (sd) and correlations (rho)
    template<class Type>
    matrix<Type> TAU( vector<Type> &sd, vector<Type> &rho ) {

        // define matrix
        int p = sd.size();
        matrix<Type> mat(p,p);
        
        // fill matrix
        int k = 0;        
        for(int i=0; i<p; i++) {
            for(int j=0; j<p; j++) { 
                if(i==j) { 
                    mat(i,i)=sd(i)*sd(i);
                } //exp(2*sdP(i));
                if(i>j){ 
                    mat(i,j)=sd(i)*sd(j)*rho(k); 
                    mat(j,i)=mat(i,j);
                    k++;
                }
            }
        }
        return mat;
    }

    // covariance matrix for relationship effects - univariate case:
    template<class Type>
    matrix<Type> uTAU_D( vector<Type> &sd, vector<Type> &rho  ) {

        // define matrix
        matrix<Type> mat(2,2);
        
        // fill matrix:
        mat(0,0)=sd(0)*sd(0);//exp(2*sdD(0));
        mat(1,1)=mat(0,0);
        mat(0,1)=sd(0)*sd(0)*rho(0);
        mat(1,0)=mat(0,1);

        return mat;

    }

    // covariance matrix for relationship effects - univariate case:
    template<class Type>
    matrix<Type> bTAU_D( vector<Type> &sd, vector<Type> &rho, int &no_var ) {

        // define matrix
        int p = sd.size();
        matrix<Type> mat(2*p,2*p);

        if ( no_var == 1 ) {
            // fill matrix:
            mat(0,0)=sd(0)*sd(0);//exp(2*sdD(0));
            mat(1,1)=mat(0,0);
            mat(0,1)=sd(0)*sd(0)*rho(0);
            mat(1,0)=mat(0,1);
        } else if ( no_var == 2 ) {
            // fill matrix:
            mat(0,0)=sd(0)*sd(0);//exp(2*sdD(0));
            mat(1,1)=mat(0,0);//exp(2*sdD(0));
            mat(2,2)=sd(1)*sd(1);//exp(2*sdD(1));
            mat(3,3)=mat(2,2);//exp(2*sdD(1));
            // cov measure 1:
            mat(0,1)=sd(0)*sd(0)*rho(0);
            mat(1,0)=mat(0,1);
            // intra:
            mat(0,2)=sd(0)*sd(1)*rho(1);
            mat(2,0)=mat(0,2);
            mat(1,3)=mat(0,2);
            mat(3,1)=mat(0,2);
            // inter:
            mat(0,3)=sd(0)*sd(1)*rho(2);
            mat(3,0)=mat(0,3);
            mat(1,2)=mat(0,3);
            mat(2,1)=mat(0,3);
            // cov measure 2:
            mat(2,3)=sd(1)*sd(1)*rho(3);
            mat(3,2)=mat(2,3);
        }

        return mat;

    }
        
}