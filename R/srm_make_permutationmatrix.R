
#---- this function generates a permutation matrix 

srm_make_permutationmatrix <- function( np, nv = 2 ) 
{
  
	# Gesamtdimension: np Personen * nv Variablen * 2 Effekte (p,t)
	n_total <- np * nv * 2

	if (nv == 1) return(diag(1, n_total))
	  
	# "Person-zuerst"-Reihenfolge: für Person i, Variable v, Effekt e
	# Index = (i-1)*nv*2 + (v-1)*2 + e
	  
	# "Variable-zuerst"-Reihenfolge: für Variable v, Person i, Effekt e  
	# Index = (v-1)*np*2 + (i-1)*2 + e
	  
	perm <- integer(n_total)
	for (i in 1:np) {
	  	for (v in 1:nv) {
	      	for (e in 1:2) {
	        	from <- (i-1)*nv*2 + (v-1)*2 + e  # Person-zuerst
	        	to   <- (v-1)*np*2 + (i-1)*2 + e  # Variable-zuerst
	        	perm[to] <- from
	      	}
	    }
	}
	  
	P <- matrix(0, n_total, n_total)
	for (j in 1:n_total) P[j, perm[j]] <- 1
	return( P )
}

