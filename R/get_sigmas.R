
#---- this function computes the covariance matrices for the SRM/TSRM components:

get_sigmas <- function( parm_list=NULL )
{
	
	#- get matrices:
	SD_G  <- parm_list[["SD_G"]]
	SD_P  <- parm_list[["SD_P"]]
	SD_D  <- parm_list[["SD_D"]]
	SD_T  <- parm_list[["SD_T"]]
	RHO_G <- parm_list[["RHO_G"]] 
	RHO_P <- parm_list[["RHO_P"]]
	RHO_D <- parm_list[["RHO_D"]]
	RHO_T <- parm_list[["RHO_T"]]
    
    #- compute sigmas
    parm_list$SIGMA_G <- if ( ncol(SD_G) > 0 ) t(SD_G) %*% RHO_G %*% SD_G else parm_list[["SIGMA_G"]]
    parm_list$SIGMA_P <- t(SD_P) %*% RHO_P %*% SD_P
    parm_list$SIGMA_D <- t(SD_D) %*% RHO_D %*% SD_D
    parm_list$SIGMA_T <- if ( ncol(SD_T) > 0 ) t(SD_T) %*% RHO_T %*% SD_T else parm_list[["SIGMA_T"]]
	
	#- consider the random group:
	return( parm_list )

} 