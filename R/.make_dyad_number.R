
#----- function that generates an identifier for the dyad and the dyad type

make_dyad_number <- function( data=NULL, p_var=NULL, g_var=NULL, dyad_name=NULL, maxg=1e2 )
{
	
	#- start adding the dyad number:
	tmp_data <- data
	groups   <- unique( tmp_data[,g_var] )
	ngroups  <- length( groups)
  	tmp_data$Dyad_type <- tmp_data$Dyad <- NA

	for ( gg in 1:ngroups) {

	    #- get group data:
	    idx      <- which( tmp_data[,g_var] == groups[ gg ] )
	    data_idx <- tmp_data[idx,]

	    #- make the first dyadic identifier: ij and ji are given the same number
	    data1 <- data_idx[,p_var]
	    sm1   <- data1[,p_var[1]] < data1[,p_var[2]]
	    data_idx[,p_var[1]]  <- maxg*( maxg + ifelse(sm1, data1[,p_var[1]], data1[,p_var[2]]))
	    data_idx[,p_var[2]]  <- maxg + ifelse(sm1, data1[,p_var[2]], data1[,p_var[1]])
	    tmp_data[idx,"Dyad"] <- data_idx[,p_var[1]] + data_idx[,p_var[2]]

	    #- make a second identifier: is 1 for ij and 2 for ji
	    tmp_data[idx,"Dyad_type"] <- ifelse(sm1, 1, 2)

	}

	if ( !is.null( dyad_name ) ) {
		name_idx <- ncol( tmp_data )
		colnames( tmp_data )[c(name_idx-1,name_idx)] <- c( dyad_name, paste0(dyad_name,"_type"))
	}

	return( tmp_data )

}