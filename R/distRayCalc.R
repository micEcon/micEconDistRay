distRayCalc <- function( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE ) {
  
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
  
    
  #### calculate angles and distances of output quantities ####
  datY <- distRayCalcYVar( yNames, data, fixThetas = fixThetas )
  
    
  #### calculate the logarithmic distance ####
  result <- rep( coefList$alpha0, nrow( data ) )
  for( i in 1:nInp ) {
    result <- result + coefList$alphaVec[i] * log( data[[ xNames[i] ]] )
  }
  for( i in 1:nOut ) {
    result <- result + coefList$betaVec[i] * datY[[ paste0( "y", i ) ]]
  }
  if( nCon > 0 ) {
    for( i in 1:nCon ) {
      result <- result + coefList$deltaVec[i] * data[[ zNames[i] ]]
    }
  }
  if( nShi > 0 ) {
    for( i in 1:nShi ) {
      result <- result + coefList$deltaVec[nCon+i] * data[[ sNames[i] ]]
    }
  }
  if( form == "tl" ) {
    for( i in 1:nInp ) {
      for( j in 1:nInp ) {
        result <- result + 0.5 * coefList$alphaMat[i,j] * 
          log( data[[ xNames[i] ]] ) * log( data[[ xNames[j] ]] )
      }
    }
    for( i in 1:nOut ) {
      for( j in 1:nOut ) {
        result <- result + 0.5 * coefList$betaMat[i,j] * 
          datY[[ paste0( "y", i ) ]] * datY[[ paste0( "y", j ) ]]
      }
    }
    for( i in 1:nInp ) {
      for( j in 1:nOut ) {
        result <- result + coefList$psiMat[i,j] * 
          log( data[[ xNames[i] ]] ) * datY[[ paste0( "y", j ) ]]
      }
    }
    if( nCon > 0 ) {
      for( i in 1:nCon ) {
        for( j in 1:nCon ) {
          result <- result + 0.5 * coefList$deltaMat[i,j] * 
            data[[ zNames[i] ]] * data[[ zNames[j] ]]
        }
      }
      for( i in 1:nInp ) {
        for( j in 1:nCon ) {
          result <- result + coefList$xiMat[i,j] * 
            log( data[[ xNames[i] ]] ) * data[[ zNames[j] ]]
        }
      }
      for( i in 1:nOut ) {
        for( j in 1:nCon ) {
          result <- result + coefList$zetaMat[i,j] * 
            datY[[ paste0( "y", i ) ]] * data[[ zNames[j] ]]
        }
      }
    }
  }
  names( result ) <- row.names( data )
  
  # return the calculated logarithmic distance value
  return( result )
}