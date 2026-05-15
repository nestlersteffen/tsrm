
#---- this function computes the covariance matrix for a single group:

srm_compute_V <- function( Zg=NULL, Zp=NULL, Zd=NULL, SIGMA_G=NULL, SIGMA_P=NULL, 
	SIGMA_D=NULL, random_group=FALSE )
{
	#- compute standard V:
	Vp <- Zp%*%( SIGMA_P )%*%t(Zp)
	Vd <- Zd%*%( SIGMA_D )%*%t(Zd)
	V  <- Vp + Vd
	#- consider the random group:
	if ( random_group ) {
		V <- V + Zg%*%( SIGMA_G )%*%t(Zg)
	}
	return( V )
} 

#---- this function computes the inverse of the covariance matrix for a single group

srm_compute_invV <- function ( Zp=NULL, Zd=NULL, SIGMA_G=NULL, SIGMA_P=NULL ) 
{
	#- we compute the Woodbury identity:
	invA  <- solve( Zd%*%( SIGMA_D )%*%t(Zd) ) 
	temp0 <- t(Zp)%*%invA 
	temp1 <- solve( SIGMA_P ) + temp0%*%Zp
	temp2 <- solve( temp1 )
	invV  <- invA - t( temp0 ) %*% temp2 %*% temp0
	return( invV )
}