
#- this function creates a list that contains the design matrices and a matrix with group infos 

srm_make_datalist <- function( srm_data = NULL, names_list = NULL, 
	args_list = NULL ) 
{

	#- step 1: we first generate a list with group specific design matrices:
  	data_list_groups <- srm_make_datalist_groups( srm_data = srm_data, names_list = names_list, 
  		args_list = args_list ) 

  	#- step 2: now we combine the matrices - for version 1, ignore this step!
  	data_list <- srm_make_datalist_combine( data_list_groups = data_list_groups, 
  		no_var = names_list$no_var )
  	return( data_list )

}