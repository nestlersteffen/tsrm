
#include "dmvnorm.h"

//- thread-safe sampling from multivariate normal
Eigen::VectorXd rmvnorm_rcpp( const Eigen::VectorXd& mean, const Eigen::MatrixXd& cov ) 
{
    int n = mean.size();
    // thread_local instead of static
    thread_local static std::mt19937 rng(
        std::random_device{}() + 
        std::hash<std::thread::id>{}(std::this_thread::get_id())
    );
    std::normal_distribution<double> norm(0.0, 1.0);
    Eigen::LLT<Eigen::MatrixXd> llt(cov);
    Eigen::MatrixXd L = llt.matrixL();
    Eigen::VectorXd z(n);
    for(int i = 0; i < n; i++) {
        z(i) = norm(rng);
    }
    return mean + L * z;
}

//- function to compute log determinant
double log_determinant_rcpp( const Eigen::MatrixXd& M ) 
{
    Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
    return ldlt.vectorD().array().log().sum();
}

//- function to compute multivariate normal log density:
double dmvnorm_rcpp( const Eigen::VectorXd& x, 
	                 const Eigen::VectorXd& mu, 
	                 const Eigen::MatrixXd& sigma,
	                 bool use_log ) 
{
    int n = x.size();
    Eigen::VectorXd resid = x - mu;  
    // Cholesky decomposition for solving
    Eigen::LDLT<Eigen::MatrixXd> ldlt(sigma);
    Eigen::VectorXd solved = ldlt.solve(resid);
    // all for the mvn:
    double quad_form = resid.dot(solved);
    double log_det = ldlt.vectorD().array().log().sum();
    double out = -0.5 * (n * std::log(2.0 * M_PI) + log_det + quad_form);
    if ( !use_log ) out = std::exp( out ); 
    return out;
}

//- function to compute multivariate normal log density:
double dmvnorm_centered_rcpp( const Eigen::VectorXd& x, 
                              const Eigen::MatrixXd& sigma,
                              bool use_log ) 
{
    int n = x.size();
    // Cholesky decomposition for solving
    Eigen::LDLT<Eigen::MatrixXd> ldlt(sigma);
    Eigen::VectorXd solved = ldlt.solve(x);
    // all for the mvn:
    double quad_form = x.dot(solved);
    double log_det = ldlt.vectorD().array().log().sum();
    double out = -0.5 * (n * std::log(2.0 * M_PI) + log_det + quad_form);
    if ( !use_log ) out = std::exp( out ); 
    return out;
}