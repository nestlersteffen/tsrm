
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
            if ( ( pos[1] == 1 & pos[2] == 3 ) | ( pos[1] == 3 & pos[2] == 5 ) |  ( pos[1] == 1 & pos[2] == 5 ) ) {
                x_nn[ pos[1] + 1, pos[2] + 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] + 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
            if ( ( pos[1] == 1 & pos[2] == 6 ) | ( pos[1] == 1 & pos[2] == 4 ) |  ( pos[1] == 3 & pos[2] == 6 ) ) {
                x_nn[ pos[1] + 1, pos[2] - 1 ] <- parm[ index_nn ]
                x_nn[ pos[2] - 1, pos[1] + 1 ] <- parm[ index_nn ]
            }
        } 
        if ( type == "SD_T" ) {
            diag( x_nn )[ c(2:6) ] <- parm[ index_nn ]
        }
        if ( type == "RHO_T") {
            x_nn[ pos[2] , pos[1] ] <- parm[ index_nn ]
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
                x_nn[1,5] <- x_nn[5,1] <- x_nn[2,3] <- x_nn[3,2] <- parm[ index_nn ]
                x_nn[2,6] <- x_nn[6,2] <- x_nn[3,6] <- x_nn[6,3] <- parm[ index_nn ]
                x_nn[4,5] <- x_nn[5,4] <- parm[ index_nn ]
            }
        }
        parm_list[[ type ]] <- x_nn       
    }
    return(parm_list)
}