distRayCalcYVar <- function( yNames, data, calcOmegas = FALSE, 
  fixThetas = FALSE ) {

  #### calculate angles and distances of output quantities ####
  nOut <- length( yNames )
  datY <- data.frame( row.names = row.names( data ) )
  for( i in 1:nOut ) {
    datY[[ paste0( "dist_", i ) ]] <- 
      sqrt( rowSums( data[ , yNames[ i:nOut ], drop = FALSE ]^2 ) )
  }
  for( i in 1:( nOut - 1 ) ) {
    datY[[ paste0( "y", i ) ]] <- acos( data[[ yNames[ i ] ]] / 
        datY[[ paste0( "dist_", i ) ]] )
    if( fixThetas ) {
      datY[[ paste0( "y", i ) ]][ is.na( datY[[ paste0( "y", i ) ]] ) &
          datY[[ paste0( "dist_", i ) ]] == 0 ] <-
        acos( 1 / sqrt( nOut + 1 - i ) )
    }
  }
  datY[[ paste0( "y", nOut ) ]] <- log( datY$dist_1 )
  
  if( calcOmegas ) {
    for( j in 1:( nOut - 1 ) ) {
      for( i in 1:nOut ) {
        if( i < j ) {
          datY[[ paste0( "omega_", j, "_", i ) ]] <- 0
        } else if( i == j ) {
          datY[[ paste0( "omega_", j, "_", i ) ]] <- 
            ifelse( data[[ yNames[ j ] ]] > 0 & 
                rowSums( data[ , yNames[ (j+1):nOut ], drop = FALSE ] ) == 0,
              0, 
              ifelse( rowSums( data[ , yNames[ j:nOut, drop = FALSE ]] ) == 0, 
              NA,
              ( data[[ yNames[ i ] ]] * data[[ yNames[ j ] ]] -
                  datY[[ paste0( "dist_", j ) ]]^2 ) /
                ( datY[[ paste0( "dist_", j ) ]]^2 * 
                    datY[[ paste0( "dist_", j + 1 ) ]] ) ) )
        } else if( i > j ) {
          datY[[ paste0( "omega_", j, "_", i ) ]] <-
            ifelse( rowSums( data[ , yNames[ j:nOut ], drop = FALSE ] ) == 0, 
              NA,
              ifelse( rowSums( data[ , yNames[ (j+1):nOut ], drop = FALSE ] ) == 
                  data[[ yNames[i] ]],
            data[[ yNames[ j ] ]] /
            datY[[ paste0( "dist_", j ) ]]^2,
                  data[[ yNames[ i ] ]] * data[[ yNames[ j ] ]] /
            ( datY[[ paste0( "dist_", j ) ]]^2 * 
                datY[[ paste0( "dist_", j + 1 ) ]] ) ) )
        }
      }
    }
  }
  
  return( datY )
}