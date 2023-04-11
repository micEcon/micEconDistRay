distRayCheckNumVars <- function( nInp, nOut, nCon, nShi, form, conDummy ) {

  # check arguments nInp, nOut, nCon, and nShi
  for( n in c( "nInp", "nOut", "nCon", "nShi" ) ) {
    i <- get( n )
    if( length(i) != 1 || !is.numeric(i) ) {
      return( paste0( "argument '", n, "' must be a single numeric value" ) )
    }
    if( i != round(i) || i < 0 ) {
      return( paste0( "argument '", n, "' must be a non-negative integer" ) )
    }
  }
  if( !( nInp > 0 && nOut > 0 ) ) {
    return( "there must be at least one input and at least one output" )
  }
  
  #### check argument 'conDummy' ####
  if( !is.null( conDummy ) ) {
    if( form == "cd" ) {
      warning( "argument 'conDummy' is ignored if argument 'form' is 'cd'" )
    } else if( nCon == 0 ) {
      warning( "argument 'conDummy' is ignored if argument 'nCon' is '0'" )
    } else {
      if( !is.numeric( conDummy ) ) {
        return( "argument 'conDummy' must be a vector of numeric values" )
      }
      if( any( conDummy < 1 ) || any( conDummy > nCon ) ) {
        return( "argument 'conDummy' must be a vector of integer values",
         "between 1 and argument 'nCon'" )
      }
    }
  }
  
  return( NULL )
}
