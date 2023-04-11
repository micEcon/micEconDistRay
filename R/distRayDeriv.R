distRayDeriv <- function( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE, 
  numDeriv = FALSE, eps = 1e-6 ) {
  
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
  
    
  #### calculate angles, distances, and "Omegas" of output quantities ####
  datY <- distRayCalcYVar( yNames, data, calcOmegas = TRUE, fixThetas = TRUE )
  
    
  #### calculate the "predicted" distance ####
  datY$distPred <- exp( distRayCalc( 
    xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
    coef = coefList, data = data, form = form, conDummy = conDummy,
    fixThetas = TRUE ) )

  
  #### create an empty data.frame to which the derivatives will be added ####
  result <- data.frame( row.names = row.names( data ) )

  ### numerical finit-difference calculations
  if( numDeriv ) {
    for( i in 1:nInp ) {
      dataEps <- data
      dataEps[[ xNames[i] ]] <- data[[ xNames[i] ]] + eps
      dataEps$distPred <- exp( distRayCalc( 
        xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
        coef = coefList, data = dataEps, form = form, conDummy = conDummy ) )
      result[[ paste0( "d_x_", i ) ]] <- 
        ( dataEps$distPred - datY$distPred ) / eps
    }
    for( i in 1:nOut ) {
      dataEps <- data
      dataEps[[ yNames[i] ]] <- data[[ yNames[i] ]] + eps
      dataEps$distPred <- exp( distRayCalc( 
        xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
        coef = coefList, data = dataEps, form = form, conDummy = conDummy ) )
      result[[ paste0( "d_y_", i ) ]] <- 
        ( dataEps$distPred - datY$distPred ) / eps
    }
    if( nCon > 0 ) {
      for( i in 1:nCon ) {
        dataEps <- data
        dataEps[[ zNames[i] ]] <- data[[ zNames[i] ]] + eps
        dataEps$distPred <- exp( distRayCalc( 
          xNames = xNames, yNames = yNames, zNames = zNames, sNames = sNames, 
          coef = coefList, data = dataEps, form = form, conDummy = conDummy ) )
        result[[ paste0( "d_z_", i ) ]] <- 
          ( dataEps$distPred - datY$distPred ) / eps
      }
    }
    
  } else {
    #### calculate the derivatives wrt the input quantities ####
    for( i in 1:nInp ) {
      result[[ paste0( "d_x_", i ) ]] <- coefList$alphaVec[ i ]
      if( form == "tl" ) {
        for( j in 1:nInp ) {
          result[[ paste0( "d_x_", i ) ]] <- result[[ paste0( "d_x_", i ) ]] +
            coefList$alphaMat[ i, j ] * log( data[[ xNames[j] ]] )
        }
        for( j in 1:nOut ) {
          result[[ paste0( "d_x_", i ) ]] <- result[[ paste0( "d_x_", i ) ]] +
            coefList$psiMat[ i, j ] * datY[[ paste0( "y", j ) ]]
        }
        if( nCon > 0 ) {
          for( j in 1:nCon ) {
            result[[ paste0( "d_x_", i ) ]] <- result[[ paste0( "d_x_", i ) ]] +
              coefList$xiMat[ i, j ] * data[[ zNames[j] ]]
          }
        }
      }
      result[[ paste0( "d_x_", i ) ]] <- result[[ paste0( "d_x_", i ) ]] *
        datY$distPred / data[[ xNames[i] ]]
    }

    #### calculate the derivatives wrt the output quantities ####
    for( i in 1:nOut ) {
      result[[ paste0( "d_y_", i ) ]] <- 0
      for( j in 1:( nOut - 1 ) ) {
        result[[ paste0( "d_y_", i ) ]] <- 
          result[[ paste0( "d_y_", i ) ]] +
          coefList$betaVec[ j ] * datY[[ paste0( "omega_", j, "_", i ) ]]
      }
      result[[ paste0( "d_y_", i ) ]] <- 
        result[[ paste0( "d_y_", i ) ]] +
        coefList$betaVec[ nOut ] * data[[ yNames[ i ] ]] / datY$dist_1^2
      if( form == "tl" ) {
        for( j in 1:( nOut - 1 ) ) {
          for( k in 1:( nOut - 1 ) ) {
            result[[ paste0( "d_y_", i ) ]] <- 
              result[[ paste0( "d_y_", i ) ]] +
              coefList$betaMat[ j, k ] * datY[[ paste0( "y", k ) ]] *
              datY[[ paste0( "omega_", j, "_", i ) ]] 
          }
        }      
        for( j in 1:( nOut - 1 ) ) { 
          result[[ paste0( "d_y_", i ) ]] <- 
            result[[ paste0( "d_y_", i ) ]] +
            coefList$betaMat[ j, nOut ] * 
            ( datY[[ paste0( "omega_", j, "_", i ) ]] * log( datY$dist_1 ) +
                datY[[ paste0( "y", j ) ]] * data[[ yNames[ i ] ]] / 
                datY$dist_1^2 )
        }
        result[[ paste0( "d_y_", i ) ]] <- 
          result[[ paste0( "d_y_", i ) ]] +
          coefList$betaMat[ nOut, nOut ] * log( datY$dist_1 ) * 
          data[[ yNames[ i ] ]] / datY$dist_1^2 
        
        for( j in 1:nInp ) {
          for( k in 1:( nOut - 1 ) ) {
            result[[ paste0( "d_y_", i ) ]] <- 
              result[[ paste0( "d_y_", i ) ]] +
              coefList$psiMat[ j, k ] * log( data[[ xNames[ j ] ]] ) *
              datY[[ paste0( "omega_", k, "_", i ) ]]
          }
        }
        for( j in 1:nInp ) {
          result[[ paste0( "d_y_", i ) ]] <- 
            result[[ paste0( "d_y_", i ) ]] +
            coefList$psiMat[ j, nOut ] * log( data[[ xNames[ j ] ]] ) * 
            data[[ yNames[ i ] ]] / datY$dist_1^2
        }
        if( nCon > 0 ) {
          for( j in 1:( nOut - 1 ) ) {
            for( k in 1:nCon ) {
              result[[ paste0( "d_y_", i ) ]] <- 
                result[[ paste0( "d_y_", i ) ]] +
                coefList$zetaMat[ j, k ] * data[[ zNames[ k ] ]] *
                datY[[ paste0( "omega_", j, "_", i ) ]] 
            }  
          }
          for( j in 1:nCon ) {
            result[[ paste0( "d_y_", i ) ]] <- 
              result[[ paste0( "d_y_", i ) ]] +
              coefList$zetaMat[ nOut, j ] * data[[ zNames[ j ] ]] * 
              data[[ yNames[ i ] ]] / datY$dist_1^2  
          }
        }
      }
      result[[ paste0( "d_y_", i ) ]] <- 
        result[[ paste0( "d_y_", i ) ]] * datY$distPred
    }
    
    #### calculate the derivatives wrt the control and shifter variables ####
    if( nCon > 0 ) {
      for( i in 1:nCon ) {
        result[[ paste0( "d_z_", i ) ]] <- coefList$deltaVec[ i ]
        if( form == "tl" ) {
          for( j in 1:nCon ) {
            result[[ paste0( "d_z_", i ) ]] <- result[[ paste0( "d_z_", i ) ]] + 
              coefList$deltaMat[ i, j ] * data[[ zNames[ j ] ]]
          }
          for( j in 1:nInp ) {
            result[[ paste0( "d_z_", i ) ]] <- result[[ paste0( "d_z_", i ) ]] + 
              coefList$xiMat[ j, i ] * log( data[[ xNames[ j ] ]] )
          }
          for( j in 1:nOut ) {
            result[[ paste0( "d_z_", i ) ]] <- result[[ paste0( "d_z_", i ) ]] + 
              coefList$zetaMat[ j, i ] * datY[[ paste0( "y", j ) ]]
          }
        }
        result[[ paste0( "d_z_", i ) ]] <- result[[ paste0( "d_z_", i ) ]] *
          datY$distPred
      }
    }
    if( nShi > 0 ) {
      for( i in 1:nShi ) {
        result[[ paste0( "d_z_", nCon + i ) ]] <-
          coefList$deltaVec[ nCon + i ]
      }
    }
  }
  
    
  # return the data.frame with the calculated derivatives
  return( result )
}