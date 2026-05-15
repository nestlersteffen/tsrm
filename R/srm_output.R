
#---- functions to print the results 

srm_print_optinfos <- function( object = NULL, digits = 3L )
{
  #- define the texts:
  if ( object$args_list$with_reml ) {
    estimator <- "REML"
  } else {
    estimator <- "ML"
  }
  texts <- c("Social relations model fit with ",
             "Log-likelihood: ",
             "Deviance: ",
             "AIC: ",
             "Number of round-robin variables: ",
             "Number of round-robin groups: ",
             "Number of parameters: ")
  values <- c( estimator,
               round( object$ll, digits),
               round( object$deviance, digits),
               round( object$aic, digits),
               object$data_list$nv,
               nrow( object$data_list$groupinfo ),
               length( unique( object$parm_table$index ) ) )
  #- print texts
  # char.format <- paste("%0s%-", 14, "s", sep = "")
  for ( i in 1:length( texts ) ) {
    # tmp <- paste( sprintf( char.format, "", texts[i] ), values[i], sep = "" )                       
    tmp <- paste0( texts[i], values[i] ) 
    cat(tmp)
    cat("\n")
  }
  cat("\n")
}

#- this function is inspired by lavaans print function

srm_print_parameters <- function( object = NULL, digits = 3L )
{
  
    #- get parmtable from object:
    parm_table <- object$parm_table
      
    #- we round the values in the table:
    parm_table[,c("est","se","z","p")] <- round( parm_table[,c("est","se","z","p")], digits )

    #- get number of round robin variables:
    nv <- object$data_list$nv

    #- print info for groups
    if ( any( parm_table$type == "SD_G" ) ) {
        group_info <- paste0("Random group parameters not shown in output.", "\n")
        group_info <- paste0( group_info, "See parm_table slot for information on these parameters" )
    }
  
    #- ------------------------------------------
    #-    Standard deviations and correlations
    #- ------------------------------------------
  
    #- person effects:
    SD_P  <- round( diag( object$parm_list$SD_P ), digits )
    VAR_P <- round( diag( object$parm_list$SD_P )^2, digits )
    RHO_P <- matrix("",nrow=nv*2,ncol=nv*2-1)
    for(i in 2:nrow(RHO_P)) {
        for(j in 1:(i-1)) {
          RHO_P[i, j] <- format( round( object$parm_list$RHO_P[i,j], digits ), nsmall = digits)
        }
    }
    if ( nv == 1 ) {
        tmpNames <- c( "Perc.", "Targ." )
    } else {
        tmpNames <- paste0( rep( c( "Perc.", "Targ."), nv ), rep(1:2, each = nv ) ) 
    }
    tmpTab <- data.frame( Name = tmpNames, Variance = format( VAR_P, nsmall = digits), 
        Std.Dev. = format( SD_P, nsmall = digits), Corr = RHO_P, stringsAsFactors = FALSE )
    colnames(tmpTab) <- c("Name", "Variance", "Std.Dev.", "Corr.", rep("", ncol(RHO_P)-1))
    cat( paste0( "Random effects person-level: " , "\n" ) )
    print( tmpTab, row.names = FALSE )
    cat("\n")

    #- relationship effects:
    SD_D  <- round( diag( object$parm_list$SD_D ), digits )
    VAR_D <- round( diag( object$parm_list$SD_D )^2, digits )
    RHO_D <- matrix("",nrow=nv*2,ncol=nv*2-1)
    for(i in 2:nrow(RHO_D)) {
        for(j in 1:(i-1)) {
          RHO_D[i, j] <- format( round( object$parm_list$RHO_D[i,j], digits ), nsmall = digits)
        }
    }
    if ( nv == 1 ) {
        tmpNames <- c( "D_ij", "D_ji" )
    } else {
        tmpNames <- paste0( rep( c( "D_ij.", "D_ji." ), nv ), rep(1:2, each = nv ) ) 
    }
    tmpTab <- data.frame( Name = tmpNames, Variance = format( VAR_D, nsmall = digits), 
        Std.Dev. = format( SD_D, nsmall = digits), Corr = RHO_D, stringsAsFactors = FALSE )
    colnames(tmpTab) <- c("Name", "Variance", "Std.Dev.", "Corr.", rep("", ncol(RHO_D)-1))
    cat( paste0( "Random effects dyad-level: " , "\n" ) )
    print( tmpTab, row.names = FALSE )
    cat("\n")

    #- ----------------------------
    #-   Location-fixed effects
    #- ----------------------------
  
    idx <- which( parm_table$type == "BETA" )
    tmpNames <- do.call("c",object$names_list$mu_preds[1])
    tmpTab <- data.frame( Name = tmpNames, 
            Value = format( parm_table$est[idx], nsmall = digits),
            Std.Error = format( parm_table$se[idx], nsmall = digits),
            z_value = format( parm_table$z[idx], nsmall = digits),
            p_value = format( parm_table$p[idx], nsmall = digits) )
    if ( nv != 1 ) {
        # no. of predictore per outcome:
        tmpNo   <- do.call("c",lapply( object$names_list$mu_preds, length ))
        # now we generate another column with the variable name
        tmpName <- as.character() 
        for ( i in seq( nv ) ) {
            tmpName <- c( tmpName, c( object$names_y[i], rep(" ", tmpNo[i] - 1 ) ) )
        }
        # now we add that column to the table
        tmpTab <- cbind( data.frame( Outcome = tmpName ), tmpTab )
        colnames(tmpTab) <- c("Outcome", "", "Value", "Std.Error", "z-value", "p-value")
    } else {
        colnames(tmpTab) <- c("", "Value", "Std.Error", "z-value", "p-value")
    }
    cat( paste0( "Fixed effects: ", "\n" ) )
    print( tmpTab, row.names = FALSE )
    cat("\n")
  
}

srm_out <- function( object = NULL, digits = 3L )
{
    srm_print_optinfos( object = object, digits = digits )
    srm_print_parameters( object = object, digits = digits )
}