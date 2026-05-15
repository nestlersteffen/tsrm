
#---- this function computes the covariance matrix for a single group:

srm_get_sigmas <- function( parm_list = NULL )
{
	#- get matrices:
	SD_G  <- parm_list[["SD_G"]]
	SD_P  <- parm_list[["SD_P"]]
	SD_D  <- parm_list[["SD_D"]]
	RHO_G <- parm_list[["RHO_G"]] 
	RHO_P <- parm_list[["RHO_P"]] 
    RHO_D <- parm_list[["RHO_D"]]  
    #- compute sigmas
    parm_list$SIGMA_G <- t(SD_G) %*% RHO_G %*% SD_G
    parm_list$SIGMA_P <- t(SD_P) %*% RHO_P %*% SD_P
    parm_list$SIGMA_D <- t(SD_D) %*% RHO_D %*% SD_D
	#- consider the random group:
	return( parm_list )
} 