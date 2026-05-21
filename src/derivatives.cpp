
#include "derivatives.h"

Eigen::MatrixXd srm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T,
    const int& type, const Eigen::VectorXi& pos )
{
    
    //- make the result matrix:
    Eigen::MatrixXd res;

    //- now, we compute the specific derivatives
    if ( type == 1 ) { 
        int p = SD_P.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        res = ONE_ij.transpose()*RHO_P*SD_P + SD_P.transpose()*RHO_P*ONE_ij;            
    }

    if ( type == 2 ) {
        int p = RHO_P.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        res = SD_P.transpose()*ONE_ij*SD_P;       
    }

    if ( type == 3 ) { 
        int p = SD_D.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(0)-1+1, pos(1)-1+1 ) = 1;
        res = ONE_ij.transpose()*RHO_D*SD_D + SD_D.transpose()*RHO_D*ONE_ij;            
    }
    
    if ( type == 4 ) {  
        int p = RHO_D.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        if ( pos(0) == 1 && pos(1) == 3 ) {
            ONE_ij( pos(0)+1, pos(1)+1 ) = 1;
            ONE_ij( pos(1)+1, pos(0)+1 ) = 1;
        }
        if ( pos(0) == 1 && pos(1) == 4 ) {
            ONE_ij( pos(0)+1, pos(1)-1 ) = 1;
            ONE_ij( pos(1)-1, pos(0)+1 ) = 1;
        }
        res = SD_D.transpose()*ONE_ij*SD_D;         
    }

    if ( type == 5 ) { 
        int p = SD_G.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        res = ONE_ij.transpose()*RHO_G*SD_G + SD_G.transpose()*RHO_G*ONE_ij;            
    }

    if ( type == 6 ) {
        int p = RHO_G.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        res = SD_G.transpose()*ONE_ij*SD_G;       
    }

    //- output
    return res;
}

Eigen::MatrixXd tsrm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& SD_T,
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
    const Eigen::MatrixXd& RHO_T, 
    const int& type, const Eigen::VectorXi& pos )
{
    
    //- make the result matrix:
    Eigen::MatrixXd res;

    //- now, we compute the specific derivatives

    if ( type == 1 ) { 
        int p = SD_P.rows();
        Eigen::MatrixXd ONE_ij = MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        res = ONE_ij.transpose()*RHO_P*SD_P + SD_P.transpose()*RHO_P*ONE_ij;            
    }

    if ( type == 2 ) {
        int p = RHO_P.rows();
        Eigen::MatrixXd ONE_ij = MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        res = SD_P.transpose()*ONE_ij*SD_P;       
    }

    if ( type == 3 ) { 
        int p = SD_D.rows();
        Eigen::MatrixXd ONE_ij = MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(0)-1+1, pos(1)-1+1 ) = 1;
        res = ONE_ij.transpose()*RHO_D*SD_D + SD_D.transpose()*RHO_D*ONE_ij;            
    }
    
    if ( type == 4 ) {  
        int p = RHO_D.rows();
        Eigen::MatrixXd ONE_ij = MatrixXd::Zero( p, p );
        ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
        ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        if ( ( pos(0) == 1 & pos(1) == 2 ) | ( pos(0) == 3 & pos(1) == 4 ) | ( pos(0) == 5 & pos(1) == 6 ) ) {
            ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
            ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
        }
        // covariance terms different effects
        if ( ( pos(0) == 1 & pos(1) == 3 ) |  ( pos(0) == 3 & pos(1) == 5 ) |  ( pos(0) == 1 & pos(1) == 5 ) ) {
            ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
            ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
            ONE_ij( pos(0)-1+1, pos(1)-1+1 ) = 1;
            ONE_ij( pos(1)-1+1, pos(0)-1+1 ) = 1;
        }
        if ( ( pos(0) == 1 & pos(1) == 6 ) | ( pos(0) == 1 & pos(1) == 4 ) | ( pos(0) == 3 & pos(1) == 6 ) ) {
            ONE_ij( pos(0)-1, pos(1)-1 ) = 1;
            ONE_ij( pos(1)-1, pos(0)-1 ) = 1;
            ONE_ij( pos(0)-1+1, pos(1)-1-1 ) = 1;
            ONE_ij( pos(1)-1-1, pos(0)-1+1 ) = 1;
        }
        res = SD_D.transpose()*ONE_ij*SD_D;         
    }

    if ( type == 7 ) { 
        int p = SD_T.rows();
        Eigen::MatrixXd ONE_ij = Eigen::MatrixXd::Identity(p, p);
        res = RHO_T*SD_T + SD_T.transpose()*RHO_T;            
    }
    
    if ( type == 8 ) {  
        int p = RHO_T.rows();
        Eigen::MatrixXd ONE_ij = MatrixXd::Zero( p, p );
        // covariance term within judge
        if ( pos(0) == 1 & pos(1) == 2 ) {
            ONE_ij(0,1) = 1; ONE_ij(1,0) = 1; ONE_ij(2,3) = 1;
            ONE_ij(3,2) = 1; ONE_ij(4,5) = 1; ONE_ij(5,4) = 1;
        }
        // covariance term within partners
        if ( pos(0) == 1 & pos(1) == 6 ) {
            ONE_ij(0,5) = 1; ONE_ij(5,0) = 1; ONE_ij(1,3) = 1; 
            ONE_ij(3,1) = 1; ONE_ij(2,4) = 1; ONE_ij(4,2) = 1;
        }
        // covariance term within actors
        if ( pos(0) == 1 & pos(1) == 3 ) {
            ONE_ij(0,2) = 1; ONE_ij(2,0) = 1; ONE_ij(1,4) = 1;
            ONE_ij(4,1) = 1; ONE_ij(3,5) = 1; ONE_ij(5,3) = 1;
        }
        // triadic terms:
        if ( pos(0) == 1 & pos(1) == 4 ) {
            ONE_ij(0,3) = 1; ONE_ij(3,0) = 1; ONE_ij(0,4) = 1;
            ONE_ij(4,0) = 1; ONE_ij(1,2) = 1; ONE_ij(2,1) = 1;
            ONE_ij(1,5) = 1; ONE_ij(5,1) = 1; ONE_ij(2,5) = 1;
            ONE_ij(5,2) = 1; ONE_ij(3,4) = 1; ONE_ij(4,3) = 1;
        }
        res = SD_T.transpose()*ONE_ij*SD_T;                
    }

    //- output
    return res;
}