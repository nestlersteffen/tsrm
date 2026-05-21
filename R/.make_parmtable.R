
#---- this function generates the parm_table that we use in optimization

make_parmtable <- function( data_list=NULL, names_list=NULL, args_list=NULL, 
    model=c("srm","tsrm") )
{
	
    #- get no. of variables:
    no_var <- names_list$no_var

    #- everthing for BETA:
    mu_preds <- names_list[["mu_preds"]]
    if ( args_list$fixed_group ) {
        # get number of groups:
        ng <- nrow( data_list$groupinfo )
        # number of predictors ( minus 1 because mu_preds contains the intercept )
        lgt_mu_preds <- ng + length( do.call( "c", mu_preds ) ) - 1
    } else {
        lgt_mu_preds <- length( do.call( "c", mu_preds ) )
    }
    
    parm_table_BETA <- data.frame( 
        type =rep("BETA", lgt_mu_preds ), 
        pos1 =seq(1,lgt_mu_preds,1), 
        pos2 =rep(1,lgt_mu_preds), 
        ntype=rep(0,lgt_mu_preds) 
    )

    #- everything for the standard deviations and correlations
	
    #- --- SRM part ---

    if ( model == "srm" ) {

        if ( no_var == 1L ) {

            parm_table_SDCOR <- data.frame( 
                type = c( rep( "SD_P", 2), rep( "RHO_P", 1),
                          rep( "SD_D", 1), rep( "RHO_D", 1) ), 
                pos1 = c(1,2, 1, 1, 1),
                pos2 = c(1,2, 2, 1, 2),
                ntype = c(rep(1,2), rep(2,1), rep(3,1), rep(4,1) )
            ) 

            if ( args_list$random_group ) {

                parm_table_GROUP <- data.frame( 
                    type=c("SD_G"), pos1=c(1), pos2=c(1), ntype=c(5) )
                parm_table_SDCOR <- rbind( parm_table_SDCOR, parm_table_GROUP )

            }

        } else if ( no_var == 2L ) {

            parm_table_SDCOR <- data.frame( 
                type = c( rep( "SD_P", 4), rep( "RHO_P", 6),
                          rep( "SD_D", 2), rep( "RHO_D", 4) ), 
                pos1 = c(1,2,3,4, 1,1,1,2,2,3, 1,3, 1,1,1,3),
                pos2 = c(1,2,3,4, 2,3,4,3,4,4, 1,3, 2,3,4,4),
                ntype = c(rep(1,4), rep(2,6), rep(3,2), rep(4,4) )
            ) 

            if ( args_list$random_group ) {

                parm_table_GROUP <- data.frame( 
                    type=c("SD_G","SD_G","RHO_G"), pos1=c(1,2,1), pos2=c(1,2,2), ntype=c(5,5,6) )
                parm_table_SDCOR <- rbind( parm_table_SDCOR, parm_table_GROUP )

            }

        }

    }

    #- --- TSRM part ---

    if ( model == "tsrm" ) {

        if ( no_var == 1L ) {

            parm_table_SDCOR <- data.frame( 
                type = c( rep( "SD_P", 3), rep( "RHO_P", 3),
                          rep( "SD_D", 3), rep( "RHO_D", 9),   
                          rep( "SD_T", 1), rep( "RHO_T", 4) ), 
                pos1 = c(1,2,3, 1,1,2,
                         1,3,5, 1,3,5,1,1,3,1,1,3,
                         1, 1,1,1,1),
                pos2 = c(1,2,3, 2,3,3,
                         1,3,5, 2,4,6,3,6,5,4,5,6,
                         1, 2,6,3,4),
                ntype = c(rep(1,3), rep(2,3), rep(3,3), rep(4,9), 7, rep(8,4) )
            )    

        } else if ( no_var == 2L ) {

            parm_table_SDCOR <- data.frame( 
                type = c( rep( "SD_P", 6), rep( "RHO_P", 6),
                          rep( "SD_D", 6), rep( "RHO_D", 18),   
                          rep( "SD_T", 2), rep( "RHO_T", 8) ), 
                pos1 = c(1,2,3,4,5,6, 1,1,2,4,4,5,
                         1,3,5,7,9,11, 1,3,5,1,1,3,1,1,3,7,9,11,7,7,9,7,7,9,
                         1,7, 1,1,1,1,7,7,7,7),
                pos2 = c(1,2,3,4,5,6, 2,3,3,5,6,6,
                         1,3,5,7,9,11, 2,4,6,3,6,5,4,5,6,8,10,12,9,12,11,10,11,12,
                         1,7, 2,6,3,4,8,12,9,10),
                ntype = c(rep(1,6), rep(2,6), rep(3,6), rep(4,18), 7,7, rep(8,8) )
            ) 

        }

    }

    #- combine parm_table_BETA with parm_table_SDCOR:
    parm_table <- rbind( parm_table_BETA, parm_table_SDCOR )

    #- some final things:
    parm_table$index <- seq(1, nrow(parm_table), 1 )
    parm_table$type  <- as.character( parm_table$type )
    	
	return( parm_table )
} 