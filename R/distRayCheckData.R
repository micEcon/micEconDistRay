distRayCheckData <- function( xNames, yNames, zNames, sNames, data ) {
  
  #### check argument 'data' ####
  #### and give an informative error message if something is incorrect
  if( !inherits( data, "data.frame" ) ) {
    return( "argument 'data' must be a data.frame" )
  }
  allNames <- c( xNames, yNames, zNames, sNames )
  missingNames <- allNames[ !allNames %in% names( data ) ]
  if( length( missingNames ) > 0 ) {
    return( paste( "the data.frame defined by argument 'data' does not include",
      "the following variables:",
      paste( missingNames, collapse = ", " ) ) )
  }

  return( NULL )
}
