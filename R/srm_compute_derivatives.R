
srm_sigma_derivatives <- function( parm_list = NULL, type = NULL, pos = NULL ) 
{
    #- some information:
    res   <- NULL
    SD_G  <- parm_list$SD_G
    SD_P  <- parm_list$SD_P
    SD_D  <- parm_list$SD_D
    RHO_G <- parm_list$RHO_G
    RHO_P <- parm_list$RHO_P
    RHO_D <- parm_list$RHO_D
  
    #- compute relevant matrices 
    if ( type == "SD_G" ) { 
        p <- dim( SD_G )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1#exp( SD_G[pos[1],pos[2]] )
        #SD_G  <- exp( SD_G )   
        res <- t( ONE_ij ) %*% SD_G + t( SD_G ) %*% ONE_ij               
    }

    if ( type == "RHO_G" ) {
        p <- dim( RHO_G )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[2],pos[1]] <- 1
        res <- t( SD_G ) %*% ONE_ij %*% SD_G           
    }

    if ( type == "SD_P" ) { 
        p <- dim( SD_P )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1#exp( SD_P[pos[1],pos[2]] )
        #diag( SD_P ) <- exp( diag( SD_P ) )   
        res <- t( ONE_ij ) %*% RHO_P %*% SD_P + t( SD_P ) %*% RHO_P %*% ONE_ij            
    }

    if ( type == "RHO_P" ) {
        #diag( SD_P ) <- exp( diag( SD_P ) )   
        p <- dim( RHO_P )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[2],pos[1]] <- 1
        res <- t( SD_P ) %*% ONE_ij %*% SD_P           
    }
      
    if ( type == "SD_D" ) { 
        p <- dim( SD_D )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        #diag( ONE_ij ) <- 1#exp( diag( SD_D ) )
        #diag( SD_D )   <- exp( diag( SD_D ) )  
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[1]+1,pos[2]+1] <- 1 
        res <- t( ONE_ij ) %*% RHO_D %*% SD_D + t( SD_D ) %*% RHO_D %*% ONE_ij            
    }

    if ( type == "RHO_D" ) {  
        # diag( SD_D ) <- exp( diag( SD_D ) ) 
        p <- dim( RHO_D )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[2],pos[1]] <- 1
        if ( pos[1] == 1 & pos[2] == 3 ) {
            ONE_ij[pos[1]+1,pos[2]+1] <- 1
            ONE_ij[pos[2]+1,pos[1]+1] <- 1
        }
        if ( pos[1] == 1 & pos[2] == 4 ) {
            ONE_ij[pos[1]+1,pos[2]-1] <- 1
            ONE_ij[pos[2]-1,pos[1]+1] <- 1
        }
        res <- t( SD_D ) %*% ONE_ij %*% SD_D           
    }

    #- output
    return( res )
}