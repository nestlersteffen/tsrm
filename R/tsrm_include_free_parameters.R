
#---- this function inserts parms in elements of parm_list

tsrm_include_free_parameters <- function( parm=NULL, parm_list=NULL, parm_table=NULL )
{
    NP <- max(parm_table$index)
    NOP <- nrow(parm_table) 
    for (nn in 1:NOP) {
        free_nn  <- parm_table[nn,]
        type     <- free_nn$type
        pos      <- c( free_nn$pos1, free_nn$pos2 )
        index_nn <- free_nn$index
        x_nn     <- parm_list[[ type ]]   
        x_nn[ pos[1] , pos[2] ] <- parm[ index_nn ]
        if ( type == "RHO_P" ) { 
            x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ] 
        }
        if ( type == "SD_D" ) {
            x_nn[ pos[1] + 1, pos[2] + 1 ] <- parm[ index_nn ]
        }
        if ( type == "RHO_D" ) {
            x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ]
            #- first block
            if ( ( pos[1] == 1 & pos[2] == 3 ) | ( pos[1] == 3 & pos[2] == 5 ) |  ( pos[1] == 1 & pos[2] == 5 ) ) {
                x_nn[ pos[1] + 1, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 6 ) | ( pos[1] == 1 & pos[2] == 4 ) |  ( pos[1] == 3 & pos[2] == 6 ) ) {
                x_nn[ pos[1] + 1, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
            #- second block
            if ( ( pos[1] == 7 & pos[2] == 9 ) | ( pos[1] == 9 & pos[2] == 11 ) |  ( pos[1] == 7 & pos[2] == 11 ) ) {
                x_nn[ pos[1] + 1, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 7 & pos[2] == 12 ) | ( pos[1] == 7 & pos[2] == 10 ) |  ( pos[1] == 9 & pos[2] == 12 ) ) {
                x_nn[ pos[1] + 1, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
            #- cross
            if ( pos[1] == 1 & ( pos[2] %in% c(7,9,11) ) ) {
                x_nn[ 2, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, 2 ] <- parm[ index_nn ]
            }
            if ( pos[1] == 1 & ( pos[2] %in% c(8,10,12) ) ) {
                x_nn[ 2, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, 2 ] <- parm[ index_nn ]
            }
            if ( pos[1] == 3 & ( pos[2] %in% c(7,9,11) ) ) {
                x_nn[ 4, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, 4 ] <- parm[ index_nn ]
            }
            if ( pos[1] == 3 & ( pos[2] %in% c(8,10,12) ) ) {
                x_nn[ 4, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, 4 ] <- parm[ index_nn ]
            }
            if ( pos[1] == 5 & ( pos[2] %in% c(7,9,11) ) ) {
                x_nn[ 6, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, 6 ] <- parm[ index_nn ]
            }
            if ( pos[1] == 5 & ( pos[2] %in% c(8,10,12) ) ) {
                x_nn[ 6, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, 6 ] <- parm[ index_nn ]
            }
        } 
        if ( type == "SD_T" ) {
            diag( x_nn )[ ( pos[1]+1):(pos[1]+5) ] <- parm[ index_nn ]
        }
        if ( type == "RHO_T") {
            x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ]
            #- first block
            if ( ( pos[1] == 1 & pos[2] == 2 ) ) {
                x_nn[3,4] <- x_nn[4,3] <- x_nn[5,6] <- x_nn[6,5] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 6 ) ) {
                x_nn[2,4] <- x_nn[4,2] <- x_nn[3,5] <- x_nn[5,3] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 3 ) ) {
                x_nn[2,5] <- x_nn[5,2] <- x_nn[4,6] <- x_nn[6,4] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 4 ) ) {
                x_nn[1,5] <- x_nn[5,1] <- parm[ index_nn ]
                x_nn[2,3] <- x_nn[3,2] <- x_nn[2,6] <- x_nn[6,2] <- parm[ index_nn ]
                x_nn[3,6] <- x_nn[6,3] <- x_nn[4,5] <- x_nn[5,4] <- parm[ index_nn ]
            }
            #- second block
            if ( ( pos[1] == 7 & pos[2] == 8 ) ) {
                x_nn[9,10] <- x_nn[10,9] <- x_nn[11,12] <- x_nn[12,11] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 7 & pos[2] == 12 ) ) {
                x_nn[8,10] <- x_nn[10,8] <- x_nn[9,11] <- x_nn[11,9] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 7 & pos[2] == 9 ) ) {
                x_nn[8,11] <- x_nn[11,8] <- x_nn[10,12] <- x_nn[12,10] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 7 & pos[2] == 10 ) ) {
                x_nn[7,11] <- x_nn[11,7] <- parm[ index_nn ]
                x_nn[8,9]  <- x_nn[9,8]  <- x_nn[8,12]  <- x_nn[12,8]  <- parm[ index_nn ]
                x_nn[9,12] <- x_nn[12,9] <- x_nn[10,11] <- x_nn[11,10] <- parm[ index_nn ]
            }
            #- cross-covariance
            if ( ( pos[1] == 1 & pos[2] == 7 ) ) { # 28
                x_nn[2,8]  <- x_nn[8,2]  <- parm[ index_nn ]
                x_nn[3,9]  <- x_nn[9,3]  <- parm[ index_nn ]
                x_nn[4,10] <- x_nn[10,4] <- parm[ index_nn ]
                x_nn[5,11] <- x_nn[11,5] <- parm[ index_nn ]
                x_nn[6,12] <- x_nn[12,6] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 9 ) ) { # 29
                x_nn[2,11] <- x_nn[11,2] <- parm[ index_nn ]
                x_nn[3,7]  <- x_nn[7,3]  <- parm[ index_nn ]
                x_nn[4,12] <- x_nn[12,4] <- parm[ index_nn ]
                x_nn[5,8]  <- x_nn[8,5]  <- parm[ index_nn ]
                x_nn[6,10] <- x_nn[10,6] <- parm[ index_nn ]
            } 
            if ( ( pos[1] == 1 & pos[2] == 12 ) ) { # 30
                x_nn[2,10] <- x_nn[10,2] <- parm[ index_nn ]
                x_nn[3,11] <- x_nn[11,3] <- parm[ index_nn ]
                x_nn[4,8]  <- x_nn[8,4]  <- parm[ index_nn ]
                x_nn[5,9]  <- x_nn[9,5]  <- parm[ index_nn ]
                x_nn[6,7]  <- x_nn[7,6]  <- parm[ index_nn ]
            } 
            if ( ( pos[1] == 1 & pos[2] == 8 ) ) { # 31
                x_nn[2,7]  <- x_nn[7,2]  <- parm[ index_nn ]
                x_nn[3,10] <- x_nn[10,3] <- parm[ index_nn ]
                x_nn[4,9]  <- x_nn[9,4]  <- parm[ index_nn ]
                x_nn[5,12] <- x_nn[12,5] <- parm[ index_nn ]
                x_nn[6,11] <- x_nn[11,6] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 10 ) ) { # 32
                x_nn[2,12] <- x_nn[12,2] <- parm[ index_nn ]
                # x_nn[2,9]  <- x_nn[9,2]  <- parm[ index_nn ]
                x_nn[4,11] <- x_nn[11,4] <- parm[ index_nn ]
                # x_nn[3,12] <- x_nn[12,3] <- parm[ index_nn ]
                x_nn[4,7]  <- x_nn[7,4]  <- parm[ index_nn ]
                x_nn[5,10] <- x_nn[10,5] <- parm[ index_nn ]
                x_nn[6,8]  <- x_nn[8,6]  <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 11 ) ) { # 33
                #x_nn[2,12] <- x_nn[12,2] <- parm[ index_nn ]
                x_nn[2,9]  <- x_nn[9,2]  <- parm[ index_nn ]
                x_nn[3,8]  <- x_nn[8,3]  <- parm[ index_nn ]
                #x_nn[4,11] <- x_nn[11,4] <- parm[ index_nn ]
                x_nn[3,12] <- x_nn[12,3] <- parm[ index_nn ]
                x_nn[5,7]  <- x_nn[7,5]  <- parm[ index_nn ]
                x_nn[6,9]  <- x_nn[9,6]  <- parm[ index_nn ]
            }
           
        }
        parm_list[[ type ]] <- x_nn       
    }
    return(parm_list)
}