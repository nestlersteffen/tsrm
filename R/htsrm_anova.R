
htsrm_anova_singlegroup <- function( group_data=NULL, np=NULL, names_list=NULL ) 
{

	#- Step 1: we compute the tsrm part
	tmp_group_data <- group_data[group_data[,"measure",] == 1,]
	tmp_names_list <- names_list
	tmp_names_list$no_var <- 1
	tsrm_part   <- tsrm_anova_singlegroup( group_data=tmp_group_data, np=np, 
		names_list=tmp_names_list, with_info=TRUE )
	tsrm_parm   <- tsrm_part[["gparm"]]
	tsrm_person <- tsrm_part[["person"]]
	tsrm_dyads  <- tsrm_part[["dyads"]]

	#- Step 2: we compute the srm part
	tmp_group_data <- group_data[group_data[,"measure",] == 2,]
	tmp_group_data$measure   <- 1
	tmp_group_data$Dyad_type <- tmp_group_data$ab_type
	tmp_names_list <- names_list
	tmp_names_list$no_var <- 1
	tmp_names_list$p_var  <- tmp_names_list$p_var[1:2]
	tmp_names_list$d_var  <- "ab"
	tmp_names_list$d_var_type <- "ab_type"
	srm_part  <- srm_anova_singlegroup( group_data=tmp_group_data, np=np, 
		names_list=tmp_names_list, with_info=TRUE )
	srm_parm   <- srm_part[["gparm"]]
	srm_person <- srm_part[["person"]]
	srm_dyads  <- srm_part[["dyads"]]
	
	#- Step 3: compute cross-parms
	psumsq <- t(tsrm_person)%*%srm_person
	dsumsq <- (np-2)*t(tsrm_dyads)%*%srm_dyads
	sumsq  <- c( psumsq[2,2], psumsq[2,1], psumsq[1,2], psumsq[1,1], psumsq[3,2], psumsq[3,1],
				 dsumsq[1,1], dsumsq[1,2], dsumsq[6,1], dsumsq[6,2], dsumsq[3,1], dsumsq[3,2] )

	#- Step 4: get S and obtain final estimates
	S <- htsrm_anova_build_crossS( N=np )
	cross_parm <- S%*%sumsq

	#- Step 5: rearrange in the order of the parm_table:
	gparm <- c( tsrm_parm, srm_parm, cross_parm)
	gparm <- gparm[ c(1,2,3,24,25, 
				      4,5,32,31, 6,30,29, 34,33, 26, 
		              7,8,9,27, 
		              10,11,12,28,13,14,15,16,17,18, 35,36,39,40,38,37,
		              19,20,21,22,23) ]
	return( gparm )

}