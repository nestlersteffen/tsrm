
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
        # covariance terms different effects
        if ( ( pos[1] == 7 & pos[2] == 9 ) |  ( pos[1] == 9 & pos[2] == 11 ) |  ( pos[1] == 7 & pos[2] == 11 ) ) {
            ONE_ij[pos[1],pos[2]] <- 1
            ONE_ij[pos[2],pos[1]] <- 1
            ONE_ij[pos[1]+1,pos[2]+1] <- 1 # (2,4) / (4,6) / (2,6)
            ONE_ij[pos[2]+1,pos[1]+1] <- 1 # (4,2) / (6,4) / (6,2)
        }
        if ( ( pos[1] == 7 & pos[2] == 12 ) | ( pos[1] == 7 & pos[2] == 10 ) | ( pos[1] == 9 & pos[2] ==  12 ) ) {
            ONE_ij[pos[1],pos[2]] <- 1
            ONE_ij[pos[2],pos[1]] <- 1
            ONE_ij[pos[1]+1,pos[2]-1] <- 1 # (2,5) / (2,3) / (4,5)
            ONE_ij[pos[2]-1,pos[1]+1] <- 1 # (5,2) / (3,2) / (4,5)
        }
        if ( pos[1] == 1 & ( pos[2] %in% c(7,9,11) ) ) {
            ONE_ij[ 2, pos[2] + 1 ] <- 1
            ONE_ij[ pos[2] + 1, 2 ] <- 1
        } 
        if ( pos[1] == 1 & ( pos[2] %in% c(8,10,12) ) ) {
            ONE_ij[ 2, pos[2] - 1 ] <- 1
            ONE_ij[ pos[2] - 1, 2 ] <- 1
        }
        if ( pos[1] == 3 & ( pos[2] %in% c(7,9,11) ) ) {
            ONE_ij[ 4, pos[2] + 1 ] <- 1
            ONE_ij[ pos[2] + 1, 4 ] <- 1
        }
        if ( pos[1] == 3 & ( pos[2] %in% c(8,10,12) ) ) {
            ONE_ij[ 4, pos[2] - 1 ] <- 1
            ONE_ij[ pos[2] - 1, 4 ] <- 1
        }
        if ( pos[1] == 5 & ( pos[2] %in% c(7,9,11) ) ) {
            ONE_ij[ 6, pos[2] + 1 ] <- 1
            ONE_ij[ pos[2] + 1, 6 ] <- 1
        }
        if ( pos[1] == 5 & ( pos[2] %in% c(8,10,12) ) ) {
            ONE_ij[ 6, pos[2] - 1 ] <- 1
            ONE_ij[ pos[2] - 1, 6 ] <- 1
        }
        res <- t( SD_D ) %*% ONE_ij %*% SD_D           
    }

    if ( type == "SD_T" ) { 
        p <- dim( SD_T )[1] 
        ONE_ij <- matrix(0, nrow=p, ncol=p)
        ONE_ij[ (pos[1]):(pos[1]+5), (pos[1]):(pos[1]+5) ] <- diag(1, 6)
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
        #- second block (Variable 2: +6):
        if ( pos[1] == 7 & pos[2] == 8 ) {
            ONE_ij[7,8] <- ONE_ij[8,7] <- ONE_ij[9,10] <- ONE_ij[10,9] <- ONE_ij[11,12] <- ONE_ij[12,11] <- 1
        }
        if ( pos[1] == 7 & pos[2] == 12 ) {
            ONE_ij[7,12] <- ONE_ij[12,7] <- ONE_ij[8,10] <- ONE_ij[10,8] <- ONE_ij[9,11] <- ONE_ij[11,9] <- 1
        }
        if ( pos[1] == 7 & pos[2] == 9 ) {
            ONE_ij[7,9] <- ONE_ij[9,7] <- ONE_ij[8,11] <- ONE_ij[11,8] <- ONE_ij[10,12] <- ONE_ij[12,10] <- 1
        }
        if ( pos[1] == 7 & pos[2] == 10 ) {
            ONE_ij[7,10] <- ONE_ij[10,7] <- ONE_ij[7,11] <- ONE_ij[11,7] <- ONE_ij[8,9] <- ONE_ij[9,8] <- 1
            ONE_ij[8,12] <- ONE_ij[12,8] <- ONE_ij[9,12] <- ONE_ij[12,9] <- ONE_ij[10,11] <- ONE_ij[11,10] <- 1
        }
        #- cross block:
        if ( pos[1] == 1 & pos[2] == 7 ) {
            ONE_ij[1,7]  <- ONE_ij[7,1]  <- 1
            ONE_ij[2,8]  <- ONE_ij[8,2]  <- ONE_ij[3,9]  <- ONE_ij[9,3]  <- 1
            ONE_ij[4,10] <- ONE_ij[10,4] <- ONE_ij[5,11] <- ONE_ij[11,5] <- 1
            ONE_ij[6,12] <- ONE_ij[12,6] <- 1
        }
        if ( pos[1] == 1 & pos[2] == 9 ) {
            ONE_ij[1,9]  <- ONE_ij[9,1]  <- 1
            ONE_ij[2,11] <- ONE_ij[11,2] <- ONE_ij[3,7]  <- ONE_ij[7,3]  <- 1
            ONE_ij[4,12] <- ONE_ij[12,4] <- ONE_ij[5,8]  <- ONE_ij[8,5]  <- 1
            ONE_ij[6,10] <- ONE_ij[10,6] <- 1
        }
        if ( pos[1] == 1 & pos[2] == 12 ) {
            ONE_ij[1,12] <- ONE_ij[12,1]  <- 1
            ONE_ij[2,10] <- ONE_ij[10,2] <- ONE_ij[3,11] <- ONE_ij[11,3] <- 1
            ONE_ij[4,8]  <- ONE_ij[8,4]  <- ONE_ij[5,9]  <- ONE_ij[9,5]  <- 1
            ONE_ij[6,7]  <- ONE_ij[7,6]  <- 1
        }
        if ( pos[1] == 1 & pos[2] == 8 ) {
            ONE_ij[1,8]  <- ONE_ij[8,1]  <- 1
            ONE_ij[2,7]  <- ONE_ij[7,2]  <- ONE_ij[3,10] <- ONE_ij[10,3] <- 1
            ONE_ij[4,9]  <- ONE_ij[9,4]  <- ONE_ij[5,12] <- ONE_ij[12,5] <- 1
            ONE_ij[6,11] <- ONE_ij[11,6] <- 1
        }
        if ( pos[1] == 1 & pos[2] == 10 ) {
            ONE_ij[1,10] <- ONE_ij[10,1]  <- 1
            ONE_ij[2,9]  <- ONE_ij[9,2]  <- ONE_ij[3,12] <- ONE_ij[12,3] <- 1
            ONE_ij[4,7]  <- ONE_ij[7,4]  <- ONE_ij[5,10] <- ONE_ij[10,5] <- 1
            ONE_ij[6,8]  <- ONE_ij[8,6]  <- 1
        }
        if ( pos[1] == 1 & pos[2] == 11 ) {
            ONE_ij[1,11]  <- ONE_ij[11,1]  <- 1
            ONE_ij[2,12] <- ONE_ij[12,2] <- ONE_ij[3,8]  <- ONE_ij[8,3]  <- 1
            ONE_ij[4,11] <- ONE_ij[11,4] <- ONE_ij[5,7]  <- ONE_ij[7,5]  <- 1
            ONE_ij[6,9]  <- ONE_ij[9,6]  <- 1
        }
        res <- t( SD_T ) %*% ONE_ij %*% SD_T                
    }

    #- output
    return( res )
}