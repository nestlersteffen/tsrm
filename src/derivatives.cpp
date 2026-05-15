
#include "derivatives.h"

Eigen::MatrixXd srm_sigma_derivatives_rcpp( 
    const Eigen::MatrixXd& SD_G,
    const Eigen::MatrixXd& SD_P,
    const Eigen::MatrixXd& SD_D, 
    const Eigen::MatrixXd& RHO_G, 
    const Eigen::MatrixXd& RHO_P,
    const Eigen::MatrixXd& RHO_D,
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