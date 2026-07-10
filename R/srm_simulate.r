
#--- function to simulate a data for the univariate and the biavraite srm

srm_simulate <- function( no_group = 1, no_members = 5, no_vars = 1, SIGMA_G = NULL, 
  SIGMA_P = NULL, SIGMA_D = NULL, MU = NULL, unique=FALSE )
{
    #- a check for the group-argument:
    if ( no_group == 1 ) random_group <- FALSE
    if ( !is.null( SIGMA_G) ) random_group <- TRUE else random_group <- FALSE
    if ( random_group ) SIGMA_G <- as.matrix( SIGMA_G )

    #- checks for number of variables:
    if ( no_vars == 1 & dim( SIGMA_G )[1] != 1 ) { 
        stop("No. of variables does not match dimension of SIGMA_G.")
    }
    if ( no_vars == 1 & dim( SIGMA_P )[1] != 2 ) {
        stop("No. of variables does not match dimension of SIGMA_P.")
    }
    if ( no_vars == 1 & dim( SIGMA_D )[1] != 2 ) {
        stop("No. of variables does not match dimension of SIGMA_D.")
    }
    if ( no_vars == 2 & dim( SIGMA_G )[1] != 2 ) { 
        stop("No. of variables does not match dimension of SIGMA_G.")
    }
    if ( no_vars == 2 & dim( SIGMA_P )[1] != 4 ) {
        stop("No. of variables does not match dimension of SIGMA_P.")
    }
    if ( no_vars == 2 & dim( SIGMA_D )[1] != 4 ) {
        stop("No. of variables does not match dimension of SIGMA_D.")
    }

    #- get no. of dyads:
    no_dyads <- ( ( no_members - 1 )* no_members )/2

    #- some arguments used later:
    ncol_SIGMA_P <- ncol( SIGMA_P )
    ncol_SIGMA_D <- ncol( SIGMA_D )

    #- an empty dataframe:
    result <- NULL

    #- let's go:
    for ( ng in seq( no_group ) ) {

        #- generate the group effect
        g_vec <- if ( no_vars == 1 ) rep(0,2) else rep(0,4)
        if ( random_group ) {
            g <- mvtnorm::rmvnorm( 1, rep(0, no_vars), SIGMA_G )
            g_vec <- if ( no_vars == 1 ) rep(g, 2) else g[ c(1,1,2,2) ]
        }

        # generate person-and dyad-values:
        pvals <- mvtnorm::rmvnorm( no_members, rep( 0, ncol_SIGMA_P ), SIGMA_P )
        dvals <- mvtnorm::rmvnorm( no_dyads, rep( 0, ncol_SIGMA_D ), SIGMA_D )

        #- set dyad_index to 1:
        dyad_index <- 1
        
        #- one or two rrs?
        if ( no_vars == 1 ) { idx_a <- c(1:2); idx_p <- c(2,1) } else { idx_a <- c(1:4); idx_p <- c(2,1,4,3) }

        #- we loop through actors i:
        for ( i in 1:(no_members-1) ) {

            #- get person i's perceiver and target effect
            personi <- pvals[i,]
            
            #- we loop across targets j:
            for ( j in (i+1):no_members ) {

                #- get person j's perceiver and target effect
                personj <- pvals[j,]
                
                #- compute dyadic judgment vector:
                yd <- MU + g_vec + personi[idx_a] + personj[idx_p] + dvals[dyad_index,]

                #- save results:
                tmp <- data.frame( Group = ng, 
                    Actor = if ( unique ) (ng-1)*no_members + c(i,j) else c(i,j), 
                    Partner = if ( unique ) (ng-1)*no_members + c(j,i) else c(j,i), 
                    Dyad = dyad_index )
                if ( no_vars == 1 ) { 
                    tmp$y <- c(yd[1],yd[2]) 
                } else { 
                   tmp$y1 <- c(yd[1],yd[2])
                   tmp$y2 <- c(yd[3],yd[4])
                }
                result <- rbind( result, tmp )

                #- increase dyad_index :
                dyad_index <- dyad_index + 1

            }
      
        }
  
    }

    return( result )

}