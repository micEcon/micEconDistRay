distRayEla <- function( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE, ... ) {
  
  #### check arguments xNames, yNames, zNames, and sNames ####
  #### and give an informative error message if something is incorrect
  testRes <- distRayCheckNames( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames )
  if( !is.null( testRes ) ) {
    stop( testRes )
  }

    
  #### numbers of inputs, outputs, control variables and shifter variables ####
  nInp <- length( xNames )
  nOut <- length( yNames )
  nCon <- length( zNames )
  nShi <- length( sNames )

  
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

    
  #### check argument coef and prepare coefList ####
  coefList <- distRayCheckPrepareCoefList( nInp = nInp, nOut = nOut, 
    nCon = nCon, nShi = nShi, coef = coef, form = form, 
    conDummy = conDummy )
  if( !is.list( coefList ) ) {
    stop( coefList )
  }
  
  
  #### check argument 'data' ####
  #### and give an informative error message if something is incorrect
  testRes <- distRayCheckData( xNames = xNames, yNames = yNames, 
    zNames = zNames, sNames = sNames, data = data )
  if( !is.null( testRes ) ) {
    stop( testRes )
  }

  
  #### calculate the derivatives ####
  deriv <- distRayDeriv( 
    xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
    coef = coefList, data = data, form = form, conDummy = conDummy,
    fixThetas = TRUE, ... )
  

  #### calculate the "predicted" distance ####
  deriv$distPred <- exp( distRayCalc( 
    xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
    coef = coefList, data = data, form = form, conDummy = conDummy,
    fixThetas = fixThetas ) )

  
  #### create an empty data.frame to which the derivatives will be added ####
  result <- data.frame( row.names = row.names( data ) )

  
  #### calculate the distance elasticities of the input quantities ####
  for( i in 1:nInp ) {
    result[[ paste0( "ela_x_", i ) ]] <- deriv[[ paste0( "d_x_", i ) ]] *
      data[[ xNames[i] ]] / deriv$distPred 
  }
  testVec <- rowSums( result[ , paste0( "ela_x_", 1:nInp ) ] )
  if( !isTRUE( all.equal( testVec, ifelse( is.na( testVec ), NA, 1 ),
    check.attributes = FALSE, tol = 1e-5 ) ) ) {
    warning( "the distance elasticities of the inputs do not sum up to one",
      " (for at least one observation)" )
  }
  rm( testVec )
  
  #### calculate the distance elasticities of the output quantities ####
  for( i in 1:nOut ) {
    result[[ paste0( "ela_y_", i ) ]] <- deriv[[ paste0( "d_y_", i ) ]] *
      data[[ yNames[i] ]] / deriv$distPred 
    result[[ paste0( "ela_y_", i ) ]][ 
      is.na( result[[ paste0( "ela_y_", i ) ]] ) &
        data[[ yNames[ i ] ]] == 0 ] <- 0
  }
  
  
  #### calculate the elasticity of scale ####
  result$ela_scale <- -1 / rowSums( result[ , paste0( "ela_y_", 1:nOut ) ] )

  
  #### calculate semi-elasticities of control variables
  #### calculate the distance elasticities of the output quantities ####
  if( nCon + nShi > 0 ) {
    for( i in 1:( nCon + nShi ) ) {
      result[[ paste0( "sela_z_", i ) ]] <- deriv[[ paste0( "d_z_", i ) ]] /
        deriv$distPred 
    }
  }
  
    
  # return the data.frame with the calculated elasticities
  return( result )
}
