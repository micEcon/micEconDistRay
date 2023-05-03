distRayEst <- function( xNames, yNames, zNames = NULL, sNames = NULL, 
  data, form = "tl", method = "sfa", fixThetas = FALSE, ... ) {
  
  #### check arguments xNames, yNames, zNames, and sNames ####
  #### and give an informative error message if something is incorrect
  testRes <- distRayCheckNames( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames )
  if( !is.null( testRes ) ) {
    stop( testRes )
  }

    
  #### check argument 'data' ####
  #### and give an informative error message if something is incorrect
  testRes <- distRayCheckData( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames, data = data )
  if( !is.null( testRes ) ) {
    stop( testRes )
  }
  
  
  #### check argument 'form' ####
  #### and give an informative error message if something is incorrect
  if( length( form ) != 1 || !is.character( form ) ) {
    stop( "argument 'form' must be a single character string" )
  }
  supportedForms <- c( "cd", "tl" )
  if( ! form %in% supportedForms ) {
    stop( "argument 'form' must be one of: ",
      paste( supportedForms, collapse = ", " ) )
  }
  rm( supportedForms )
  
  
  #### check argument 'method' ####
  #### and give an informative error message if something is incorrect
  if( length( method ) != 1 || !is.character( method ) ) {
    stop( "argument 'method' must be a single character string" )
  }
  supportedMethods <- c( "lm", "sfa" )
  if( ! method %in% supportedMethods ) {
    stop( "argument 'method' must be one of: ",
      paste( supportedMethods, collapse = ", " ) )
  }
  rm( supportedMethods )

  #### numbers of inputs, outputs, control variables and shifter variables ####
  nInp <- length( xNames )
  nOut <- length( yNames )
  nCon <- length( zNames )
  nShi <- length( sNames )

  #### create data set for the estimation ####
  datEst <- distRayPrepareData( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames,
    data = data, form = form, fixThetas = fixThetas )
  
  #### create formula object for the estimation ####
  # check which control variables are dummy variables
  conDummy <- NULL
  if( nCon > 0 && form == "tl" ) {
    for( i in 1:nCon ) {
      if( length( unique( data[[ zNames[i] ]] ) ) == 2 ) {
        conDummy <- c( conDummy, i )
      }
    }
  }
  formulaEst <- distRayMakeFormula( nInp = nInp, nOut = nOut, nCon = nCon, 
    nShi = nShi, form = form, conDummy = conDummy )
  
  
  #### estimate the stochastic ray distance function ####
  if( method == "lm" ) {
    result <- lm( formulaEst, data = datEst, ... )
  } else if( method == "sfa" ) {
    result <- sfacross( formulaEst, data = datEst, ... )
  }
  
  # create a list with all components of the coefficients as separate components
  result$coefList <- distRayPrepareCoef( nInp = nInp, nOut = nOut, 
    nCon = nCon, nShi = nShi, coefVec = coef( result ), form = form, 
    conDummy = conDummy ) 
  
  # calculate the distance elasticities
  result$ela <- distRayEla( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames, coef = result$coefList, 
    data = data, form = form )
    
  # check monotonicity
  deriv <- distRayDeriv( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames, coef = result$coefList, 
    data = data, form = form )
  names( deriv ) <- sub( "^d_", "", names( deriv ) )
  result$mono <- as.data.frame( cbind( 
    deriv[ , 1:nInp ] >= 0, 
    deriv[ , (nInp+1):(nInp+nOut)] <= 0,
    all_x = rowSums( deriv[ , 1:nInp ] < 0 ) == 0,
    all_y = rowSums( deriv[ , (nInp+1):(nInp+nOut)] > 0 ) == 0 ) )
  result$mono <- cbind( result$mono, 
    all = ( result$mono$all_x & result$mono$all_y ) )
    
  # return the estimation results
  return( result )
}
