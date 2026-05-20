
tsrm_sigma_derivatives <- function( parm_list = NULL, type = NULL, pos = NULL ) 
{
    #- some information:
    res   <- NULL
    SD_P  <- parm_list$SD_P
    SD_D  <- parm_list$SD_D
    SD_T  <- parm_list$SD_T
    RHO_P <- parm_list$RHO_P
    RHO_D <- parm_list$RHO_D
    RHO_T <- parm_list$RHO_T    

    if ( type == "SD_P" ) { 
        p <- dim( SD_P )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1 
        res <- t( ONE_ij ) %*% RHO_P %*% SD_P + t( SD_P ) %*% RHO_P %*% ONE_ij            
    }

    if ( type == "RHO_P" ) {
        p <- dim( RHO_P )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[2],pos[1]] <- 1
        res <- t( SD_P ) %*% ONE_ij %*% SD_P           
    }

    if ( type == "SD_D" ) { 
        p <- dim( SD_D )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[1]+1,pos[2]+1] <- 1 
        res <- t( ONE_ij ) %*% RHO_D %*% SD_D + t( SD_D ) %*% RHO_D %*% ONE_ij            
    }

    if ( type == "RHO_D" ) {  
        p <- dim( RHO_D )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[pos[1],pos[2]] <- 1
        ONE_ij[pos[2],pos[1]] <- 1
        if ( ( pos[1] == 1 & pos[2] == 2 ) | ( pos[1] == 3 & pos[2] == 4 ) | ( pos[1] == 5 & pos[2] == 6 ) ) {
            ONE_ij[pos[1],pos[2]] <- 1
            ONE_ij[pos[2],pos[1]] <- 1
        }
        # covariance terms different effects
        if ( ( pos[1] == 1 & pos[2] == 3 ) |  ( pos[1] == 3 & pos[2] == 5 ) |  ( pos[1] == 1 & pos[2] == 5 ) ) {
            ONE_ij[pos[1],pos[2]] <- 1
            ONE_ij[pos[2],pos[1]] <- 1
            ONE_ij[pos[1]+1,pos[2]+1] <- 1 # (2,4) / (4,6) / (2,6)
            ONE_ij[pos[2]+1,pos[1]+1] <- 1 # (4,2) / (6,4) / (6,2)
        }
        if ( ( pos[1] == 1 & pos[2] == 6 ) | ( pos[1] == 1 & pos[2] == 4 ) | ( pos[1] == 3 & pos[2] == 6 ) ) {
            ONE_ij[pos[1],pos[2]] <- 1
            ONE_ij[pos[2],pos[1]] <- 1
            ONE_ij[pos[1]+1,pos[2]-1] <- 1 # (2,5) / (2,3) / (4,5)
            ONE_ij[pos[2]-1,pos[1]+1] <- 1 # (5,2) / (3,2) / (4,5)
        }
        res <- t( SD_D ) %*% ONE_ij %*% SD_D           
    }

    if ( type == "SD_T" ) { 
        p <- dim( SD_T )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        diag( ONE_ij ) <- 1
        res <- t( ONE_ij ) %*% RHO_T %*% SD_T + t( SD_T ) %*% RHO_T %*% ONE_ij            
    }

    if ( type == "RHO_T" ) {  
        p <- dim( RHO_T )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        # covariance term within judge
        if ( pos[1] == 1 & pos[2] == 2 ) {
            ONE_ij[1,2] <- ONE_ij[2,1] <- ONE_ij[3,4] <- ONE_ij[4,3] <- ONE_ij[5,6] <- ONE_ij[6,5] <- 1
        }
        # covariance term within partners
        if ( pos[1] == 1 & pos[2] == 6 ) {
            ONE_ij[1,6] <- ONE_ij[6,1] <- ONE_ij[2,4] <- ONE_ij[4,2] <- ONE_ij[3,5] <- ONE_ij[5,3] <- 1
        }
        # covariance term within actors
        if ( pos[1] == 1 & pos[2] == 3 ) {
            ONE_ij[1,3] <- ONE_ij[3,1] <- ONE_ij[2,5] <- ONE_ij[5,2] <- ONE_ij[4,6] <- ONE_ij[6,4] <- 1
        }
        # triadic terms:
        if ( pos[1] == 1 & pos[2] == 4 ) {
            ONE_ij[1,4] <- ONE_ij[4,1] <- ONE_ij[1,5] <- ONE_ij[5,1] <- ONE_ij[2,3] <- ONE_ij[3,2] <- 1
            ONE_ij[2,6] <- ONE_ij[6,2] <- ONE_ij[3,6] <- ONE_ij[6,3] <- ONE_ij[4,5] <- ONE_ij[5,4] <- 1
        }
        res <- t( SD_T ) %*% ONE_ij %*% SD_T                
    }

    #- output
    return( res )
}