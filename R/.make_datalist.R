
#- this function creates a list that contains the design matrices and a matrix with group infos 

make_datalist <- function( data=NULL, names_list=NULL, args_list=NULL, model=c("srm","tsrm") ) 
{

	#- Step 1: get model
	model <- match.arg( model )

	#- Step 2: select make_datalist_groups function depending on model:
    make_datalist_groups <- switch( model,
        "srm"  = srm_make_datalist_groups,
        "tsrm" = tsrm_make_datalist_groups,
        stop( paste( "The srm or the tsrm can be estimated." ) )
    )

	#- step 1: we first generate a list with group specific design matrices:
  	data_list_groups <- make_datalist_groups( data=data, names_list=names_list, 
  		args_list=args_list ) 

  	#- step 2: now we combine the matrices - for version 1, ignore this step!
  	data_list <- make_datalist_combine( data_list_groups=data_list_groups, 
  		no_var=names_list$no_var, model=model )
  	
  	return( data_list )

}