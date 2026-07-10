
#----- function that generates an identifier for the triad and the triad type

make_triad_number <- function( data=NULL, p_var=NULL, g_var=NULL, maxg=1e2 )
{

    #- start adding the dyad number:
    tmp_data <- data
    groups   <- unique( tmp_data[,g_var] )
    ngroups  <- length( groups )
    tmp_data$Triad_type <- tmp_data$Triad <- NA
    
    for ( gg in 1:ngroups ) {
        
        #- get group data:
        idx      <- which( tmp_data[,g_var] == groups[ gg ] )
        data_idx <- tmp_data[idx,]
        
        #- für jede Zeile: erstelle Triaden-Nummer und Typ
        for ( ii in 1:nrow(data_idx) ) {
            
            #- die drei Personen in der Triade:
            judge   <- data_idx[ii, p_var[3]]  # j.id
            actor   <- data_idx[ii, p_var[1]]  # p.id
            partner <- data_idx[ii, p_var[2]]  # t.id
            
            #- sortiere die drei Personen für eindeutige Triaden-ID:
            three_persons <- sort( c(judge, actor, partner) )
            
            #- erstelle eindeutige Triaden-Nummer (wie bei Dyaden, nur mit 3 Personen):
            triad_id <- maxg^2 * (maxg^2 + three_persons[1]) + 
                        maxg * (maxg + three_persons[2]) + 
                        three_persons[3]
            
            # Bestimme den Triade-Typ (1-6):
            # Sortierung nach (Judge, Actor, Partner)
            config <- c(judge, actor, partner)
            
            # Vergleiche mit allen 6 möglichen Permutationen der sortierten Triade:
            # perms <- rbind(
            #     c(three_persons[1], three_persons[2], three_persons[3]),  # Typ 1
            #     c(three_persons[1], three_persons[3], three_persons[2]),  # Typ 2
            #     c(three_persons[2], three_persons[1], three_persons[3]),  # Typ 3
            #     c(three_persons[2], three_persons[3], three_persons[1]),  # Typ 4
            #     c(three_persons[3], three_persons[1], three_persons[2]),  # Typ 5
            #     c(three_persons[3], three_persons[2], three_persons[1])   # Typ 6
            # )
            # perms <- perms[!duplicated(perms), ]

            perms <- rbind(
                c(three_persons[3], three_persons[1], three_persons[2]),  # Typ 1 = ijk
                c(three_persons[3], three_persons[2], three_persons[1]),  # Typ 2 = jik
                c(three_persons[2], three_persons[1], three_persons[3]),  # Typ 3 = ikj
                c(three_persons[2], three_persons[3], three_persons[1]),  # Typ 4 = kij
                c(three_persons[1], three_persons[2], three_persons[3]),  # Typ 5 = jki
                c(three_persons[1], three_persons[3], three_persons[2])   # Typ 6 = kji
            )
            perms <- perms[!duplicated(perms), ]
            
            # Finde welcher Typ:
            triad_type <- which( apply(perms, 1, function(x) all(x == config)) )
            
            # Speichere:
            tmp_data[idx[ii], "Triad"]      <- triad_id
            tmp_data[idx[ii], "Triad_type"] <- triad_type
        }
    }
        
    return( tmp_data )
}
