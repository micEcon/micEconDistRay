distRayCheckNames <- function( xNames, yNames, zNames, sNames ) {
  
  #### check arguments xNames, yNames, zNames, and sNames ####
  #### and give an informative error message if something is incorrect
  if( !inherits( xNames, "character" ) ) {
    return( "argument 'xNames' must be a vector of character strings" )
  }
  if( length( xNames ) == 0 ) {
    return( "argument 'xNames' must include at least one input quantity" )
  }
  if( !inherits( yNames, "character" ) ) {
    return( "argument 'yNames' must be a vector of character strings" )
  }
  if( length( yNames ) == 0 ) {
    return( "argument 'yNames' must include at least one output quantity" )
  }
  if( !is.null( zNames ) ) {
    if( !inherits( zNames, "character" ) ) {
      return( "argument 'zNames' must be a vector of character strings" )
    }
  }
  if( !is.null( sNames ) ) {
    if( !inherits( sNames, "character" ) ) {
      return( "argument 'sNames' must be a vector of character strings" )
    }
  }
  return( NULL )
}
