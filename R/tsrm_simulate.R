
#--- function to simulate data for the univariate and the bivaraite tsrm

tsrm_simulate <- function( no_group = 1, no_members = 3, no_vars = 1,
    SIGMA_P = NULL,  # univariate - 3x3: (A, B, C)
    SIGMA_D = NULL,  # univariate - 6x6: (AB_ij, AB_ji, AC_ij, AC_ji, BC_ij, BC_ji)
    SIGMA_T = NULL,  # univariate - 6x6: (ABC_ijk, ABC_jik, ABC_ikj, ABC_kij, ABC_jki, ABC_kji)
    MU = NULL )
{
    
    #- checks for number of variables:
    if ( no_vars == 1 & dim( SIGMA_P )[1] != 3 ) {
        stop("No. of variables does not match dimension of SIGMA_P.")
    }
    if ( no_vars == 1 & dim( SIGMA_D )[1] != 6 ) {
        stop("No. of variables does not match dimension of SIGMA_D.")
    }
    if ( no_vars == 1 & dim( SIGMA_T )[1] != 6 ) { 
        stop("No. of variables does not match dimension of SIGMA_T.")
    }
    
    if ( no_vars == 2 & dim( SIGMA_P )[1] != 6 ) {
        stop("No. of variables does not match dimension of SIGMA_P.")
    }
    if ( no_vars == 2 & dim( SIGMA_D )[1] != 12 ) {
        stop("No. of variables does not match dimension of SIGMA_D.")
    }
    if ( no_vars == 2 & dim( SIGMA_T )[1] != 12 ) { 
        stop("No. of variables does not match dimension of SIGMA_T.")
    }

    #- make mean vector:
    if ( no_vars == 2 && length(MU) != 2 ) {
        stop("Length of MU does not match no. of variables.")
    }

    #- number of triads in each group:
    no_triads <- choose(no_members, 3)
    
    #- number of dyads in each group:
    no_dyads <- choose(no_members, 2)
    
    #- dimensions:
    ncol_SIGMA_P <- ncol(SIGMA_P)  
    ncol_SIGMA_D <- ncol(SIGMA_D)  
    ncol_SIGMA_T <- ncol(SIGMA_T)  
    
    #- empty result dataframe:
    result <- NULL
    
    #- loop through groups:
    for (ng in seq(no_group)) {
        
        #- generate person effects (A, B, C for each person):
        pvals <- mvtnorm::rmvnorm(no_members, rep(0, ncol_SIGMA_P), SIGMA_P)
        
        #- generate dyad effects:
        dvals <- mvtnorm::rmvnorm(no_dyads, rep(0, ncol_SIGMA_D), SIGMA_D)
        
        #- generate triad effects:
        tvals <- mvtnorm::rmvnorm(no_triads, rep(0, ncol_SIGMA_T), SIGMA_T)
        
        #- create dyad mapping:
        dyad_map <- data.frame( p1 = integer(), p2 = integer(), dyad_id = integer() )
        dyad_id  <- 1
        for (i in 1:(no_members-1)) {
            for (j in (i+1):no_members) {
                dyad_map <- rbind(dyad_map, data.frame(p1 = i, p2 = j, dyad_id = dyad_id) )
                dyad_id <- dyad_id + 1
            }
        }
        
        #- loop through triads:
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
                    
                    #- compute all six judgments in the triad for the first variable:
                    
					# 1. y_ijk: i beurteilt j, k beobachtet
					y_ijk <- MU[1] + 
					    pvals[i, 1] + pvals[j, 2] + pvals[k, 3] +  # A_i, B_j, C_k
					    dvals[dyad_ij, 1] +  # AB_ij
					    dvals[dyad_ik, 3] +  # AC_ik
					    dvals[dyad_jk, 5] +  # BC_jk
					    teff[1]              # ABC_ijk

					# 2. y_jik: j beurteilt i, k beobachtet
					y_jik <- MU[1] + 
					    pvals[j, 1] + pvals[i, 2] + pvals[k, 3] +  # A_j, B_i, C_k
					    dvals[dyad_ij, 2] +  # AB_ji
					    dvals[dyad_jk, 3] +  # AC_jk
					    dvals[dyad_ik, 5] +  # BC_ik
					    teff[2]              # ABC_jik

					# 3. y_ikj: i beurteilt k, j beobachtet
					y_ikj <- MU[1] + 
					    pvals[i, 1] + pvals[k, 2] + pvals[j, 3] +  # A_i, B_k, C_j
					    dvals[dyad_ik, 1] +  # AB_ik
					    dvals[dyad_ij, 3] +  # AC_ij
					    dvals[dyad_jk, 6] +  # BC_kj
					    teff[3]              # ABC_ikj

					# 4. y_kij: k beurteilt i, j beobachtet
					y_kij <- MU[1] + 
					    pvals[k, 1] + pvals[i, 2] + pvals[j, 3] +  # A_k, B_i, C_j
					    dvals[dyad_ik, 2] +  # AB_ki
					    dvals[dyad_jk, 4] +  # AC_kj
					    dvals[dyad_ij, 5] +  # BC_ij
					    teff[4]              # ABC_kij

					# 5. y_jki: j beurteilt k, i beobachtet
					y_jki <- MU[1] + 
					    pvals[j, 1] + pvals[k, 2] + pvals[i, 3] +  # A_j, B_k, C_i
					    dvals[dyad_jk, 1] +  # AB_jk
					    dvals[dyad_ij, 4] +  # AC_ji
					    dvals[dyad_ik, 6] +  # BC_ki
					    teff[5]              # ABC_jki

					# 6. y_kji: k beurteilt j, i beobachtet
					y_kji <- MU[1] + 
					    pvals[k, 1] + pvals[j, 2] + pvals[i, 3] +  # A_k, B_j, C_i
					    dvals[dyad_jk, 2] +  # AB_kj
					    dvals[dyad_ik, 4] +  # AC_ki
					    dvals[dyad_ij, 6] +  # BC_ji
					    teff[6]              # ABC_kji

                    #- now, we compute the second variable
                    if ( no_vars == 2 ) {

                        # 1. y_ijk: i beurteilt j, k beobachtet
                        y2_ijk <- MU[2] + 
                            pvals[i, 4] + pvals[j, 5] + pvals[k, 6] +  # A_i, B_j, C_k
                            dvals[dyad_ij, 7] +  # AB_ij
                            dvals[dyad_ik, 9] +  # AC_ik
                            dvals[dyad_jk, 11] +  # BC_jk
                            teff[7]              # ABC_ijk

                        # 2. y_jik: j beurteilt i, k beobachtet
                        y2_jik <- MU[2] + 
                            pvals[j, 4] + pvals[i, 5] + pvals[k, 6] +  # A_j, B_i, C_k
                            dvals[dyad_ij, 8] +  # AB_ji
                            dvals[dyad_jk, 9] +  # AC_jk
                            dvals[dyad_ik, 11] +  # BC_ik
                            teff[8]              # ABC_jik

                        # 3. y_ikj: i beurteilt k, j beobachtet
                        y2_ikj <- MU[2] + 
                            pvals[i, 4] + pvals[k, 5] + pvals[j, 6] +  # A_i, B_k, C_j
                            dvals[dyad_ik, 7] +  # AB_ik
                            dvals[dyad_ij, 9] +  # AC_ij
                            dvals[dyad_jk, 12] +  # BC_kj
                            teff[9]              # ABC_ikj

                        # 4. y_kij: k beurteilt i, j beobachtet
                        y2_kij <- MU[2] + 
                            pvals[k, 4] + pvals[i, 5] + pvals[j, 6] +  # A_k, B_i, C_j
                            dvals[dyad_ik, 8] +  # AB_ki
                            dvals[dyad_jk, 10] +  # AC_kj
                            dvals[dyad_ij, 11] +  # BC_ij
                            teff[10]              # ABC_kij

                        # 5. y_jki: j beurteilt k, i beobachtet
                        y2_jki <- MU[2] + 
                            pvals[j, 4] + pvals[k, 5] + pvals[i, 6] +  # A_j, B_k, C_i
                            dvals[dyad_jk, 7] +  # AB_jk
                            dvals[dyad_ij, 10] +  # AC_ji
                            dvals[dyad_ik, 12] +  # BC_ki
                            teff[11]              # ABC_jki

                        # 6. y_kji: k beurteilt j, i beobachtet
                        y2_kji <- MU[2] + 
                            pvals[k, 4] + pvals[j, 5] + pvals[i, 6] +  # A_k, B_j, C_i
                            dvals[dyad_jk, 8] +  # AB_kj
                            dvals[dyad_ik, 10] +  # AC_ki
                            dvals[dyad_ij, 12] +  # BC_ji
                            teff[12] 

                    }

                    #- save all judgments:
                    tmp <- data.frame(
                        Group = ng, 
                        Actor = c(i, j, i, k, j, k),
                        Partner = c(j, i, k, i, k, j),
                        Judge = c(k, k, j, j, i, i),
                        Triad = triad_index,
                        Dyad_AB = c(dyad_ij, dyad_ij, dyad_ik, dyad_ik, dyad_jk, dyad_jk),
                        Dyad_AC = c(dyad_ik, dyad_jk, dyad_ij, dyad_jk, dyad_ij, dyad_ik),
                        Dyad_BC = c(dyad_jk, dyad_ik, dyad_jk, dyad_ij, dyad_ik, dyad_ij) )
                    if ( no_vars == 1 ) tmp$y <- c(y_ijk, y_jik, y_ikj, y_kij, y_jki, y_kji)
                    if ( no_vars == 2 ) {
                        tmp$y1 <- c(y_ijk, y_jik, y_ikj, y_kij, y_jki, y_kji)
                        tmp$y2 <- c(y2_ijk, y2_jik, y2_ikj, y2_kij, y2_jki, y2_kji)
                    }
                    result <- rbind(result, tmp )
                    
                    triad_index <- triad_index + 1
                }
            }
        }
    }
    
    return(result)
}