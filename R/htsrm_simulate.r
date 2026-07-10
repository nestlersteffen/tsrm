
#--- function to simulate a data for the univariate and the biavraite srm

htsrm_simulate <- function( no_group = 1, no_members = 5,  
  SIGMA_P = NULL, SIGMA_D = NULL, SIGMA_T = NULL, MU = NULL )
{
    #- some checks:
    if ( dim( SIGMA_P )[1] != 5 ) {
        stop("The dimension of SIGMA_P does not match the htsrm.")
    }
    if ( dim( SIGMA_D )[1] != 8 ) {
        stop("The dimension of SIGMA_D does not match the htsrm.")
    }
    if ( dim( SIGMA_T )[1] != 6 ) {
        stop("The dimension of SIGMA_T does not match the htsrm.")
    }
    if ( length(MU) != 2 ) {
        stop("Length of MU does not match the htsrm.")
    }

    #- number of triads in each group:
    no_triads <- choose(no_members, 3)

    #- get no. of dyads:
    no_dyads <- ( ( no_members - 1 )* no_members )/2

    #- some arguments used later:
    ncol_SIGMA_P <- ncol( SIGMA_P ) # sollte 5 sein
    ncol_SIGMA_D <- ncol( SIGMA_D ) # sollte 8 sein
    ncol_SIGMA_T <- ncol( SIGMA_T ) # sollte 6 sein

    #- an empty dataframe:
    result_triad <- result_dyad <- NULL

    #- let's go:
    for ( ng in seq( no_group ) ) {

        #- generate person effects (A, B, C for each person):
        pvals <- mvtnorm::rmvnorm(no_members, rep(0, ncol_SIGMA_P), SIGMA_P)
        
        #- generate dyad effects:
        dvals <- mvtnorm::rmvnorm(no_dyads, rep(0, ncol_SIGMA_D), SIGMA_D)
        
        #- generate triad effects:
        tvals <- mvtnorm::rmvnorm(no_triads, rep(0, ncol_SIGMA_T), SIGMA_T)

        #- ----- generate the triadic judgment

        #- step 1: create the dyad mapping for the triadic judgment
        dyad_map <- data.frame( p1 = integer(), p2 = integer(), dyad_id = integer() )
        dyad_id  <- 1
        for (i in 1:(no_members-1)) {
            for (j in (i+1):no_members) {
                dyad_map <- rbind(dyad_map, data.frame(p1 = i, p2 = j, dyad_id = dyad_id) )
                dyad_id <- dyad_id + 1
            }
        }

        #- step 2: loop through triads to build the judgments
        triad_index <- 1
        for (i in 1:(no_members-2)) {
            for (j in (i+1):(no_members-1)) {
                for (k in (j+1):no_members) {
                    
                    #- get dyad IDs for this triad:
                    dyad_ij <- dyad_map$dyad_id[dyad_map$p1 == i & dyad_map$p2 == j]
                    dyad_ik <- dyad_map$dyad_id[dyad_map$p1 == i & dyad_map$p2 == k]
                    dyad_jk <- dyad_map$dyad_id[dyad_map$p1 == j & dyad_map$p2 == k]
                    
                    #- get triad effects (for six or twelve judgments ):
                    teff <- tvals[triad_index, ]
                    
                    #- compute the six triadic judgments:
                    y_ijk <- MU[1] + pvals[i, 1] + pvals[j, 2] + pvals[k, 3] +  
                             dvals[dyad_ij, 1] + dvals[dyad_ik, 3] + dvals[dyad_jk, 5] + teff[1]
                    y_jik <- MU[1] + pvals[j, 1] + pvals[i, 2] + pvals[k, 3] +  
                             dvals[dyad_ij, 2] + dvals[dyad_jk, 3] + dvals[dyad_ik, 5] + teff[2]
                    y_ikj <- MU[1] + pvals[i, 1] + pvals[k, 2] + pvals[j, 3] +  
                             dvals[dyad_ik, 1] + dvals[dyad_ij, 3] + dvals[dyad_jk, 6] + teff[3]
                    y_kij <- MU[1] + pvals[k, 1] + pvals[i, 2] + pvals[j, 3] +  
                             dvals[dyad_ik, 2] + dvals[dyad_jk, 4] + dvals[dyad_ij, 5] + teff[4]
                    y_jki <- MU[1] + pvals[j, 1] + pvals[k, 2] + pvals[i, 3] +  
                             dvals[dyad_jk, 1] + dvals[dyad_ij, 4] + dvals[dyad_ik, 6] + teff[5]
                    y_kji <- MU[1] + pvals[k, 1] + pvals[j, 2] + pvals[i, 3] + 
                             dvals[dyad_jk, 2] + dvals[dyad_ik, 4] + dvals[dyad_ij, 6] + teff[6]

                    #- make a dataframe with these judgments:
                    result_triad <- rbind( result_triad, data.frame(
                        Group = ng, 
                        Actor = c(i, j, i, k, j, k),
                        Partner = c(j, i, k, i, k, j),
                        Judge = c(k, k, j, j, i, i),
                        y1 = c(y_ijk, y_jik, y_ikj, y_kij, y_jki, y_kji),
                        y2 = NA ) )
                    
                    triad_index <- triad_index + 1
                }
            }
        }

        #- ----- generate the dyadic judgment

        dyad_index <- 1
        for ( i in 1:(no_members-1) ) {
            for ( j in (i+1):no_members ) {

                #- get dyad IDs for this dyad:
                dyad_ij <- dyad_map$dyad_id[dyad_map$p1 == i & dyad_map$p2 == j]
                
                #- compute dyadic judgment vector:
                y_ij <- MU[2] + pvals[i,4] + pvals[j,5] + dvals[dyad_ij,7]
                y_ji <- MU[2] + pvals[j,4] + pvals[i,5] + dvals[dyad_ij,8]

                #- save results:
                result_dyad <- rbind( result_dyad, data.frame(
                        Group = ng, 
                        Actor = c(i, j),
                        Partner = c(j, i),
                        Judge = c( i, j),
                        y1 = NA,
                        y2 = c( y_ij, y_ji ) ) )

                #- increase dyad_index :
                dyad_index <- dyad_index + 1

            }
      
        }
  
    }

    #- ---- merge the two data_frames:
    result <- rbind( result_triad, result_dyad )
    return( result )

}