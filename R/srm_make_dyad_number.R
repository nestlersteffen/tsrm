
#--- function to add a dyad number:

srm_make_dyad_number <- function( srm_data = NULL, p_var = NULL, g_var = NULL )
{
	
	groups  <- unique( srm_data[,g_var] )
	ngroups <- length( groups)
  srm_data$Dyad <- srm_data$Dyad_type <- NA

  maxg <- 1e2

  for ( gg in 1:ngroups) {

    #- get group data:
    idx <- which( srm_data[,g_var] == groups[ gg ] )
    data_idx <- srm_data[idx,]

    #- make the first dyadic identifier: ij and ji are given the same number
    data1 <- data_idx[,p_var]
    sm1   <- data1[,p_var[1]] < data1[,p_var[2]]
    data_idx[,p_var[1]] <- maxg*( maxg + ifelse(sm1, data1[,p_var[1]], data1[,p_var[2]]))
    data_idx[,p_var[2]] <- maxg + ifelse(sm1, data1[,p_var[2]], data1[,p_var[1]])
    srm_data[idx,"Dyad"] <- data_idx[,p_var[1]] + data_idx[,p_var[2]]

    #- make a second identifier: is 1 for ij and 2 for ji
    srm_data[idx,"Dyad_type"] <- ifelse(sm1, 1, 2)

  }

  return( srm_data )

}