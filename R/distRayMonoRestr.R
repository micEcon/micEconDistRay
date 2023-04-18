distRayMonoRestr <- function( xNames, yNames, 
  zNames = NULL, sNames = NULL, data, form = "tl", conDummy = NULL ) {
  
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
  
  
  #### numbers of inputs, outputs, control variables and shifter variables ####
  nInp <- length( xNames )
  nOut <- length( yNames )
  nCon <- length( zNames )
  nShi <- length( sNames )
  

  #### calculate angles and distances of output quantities ####
  datY <- distRayCalcYVar( yNames, data, calcOmegas = TRUE )
  

  #### obtain names of estimated coefficients ####  
  formulaEst <- distRayMakeFormula( nInp = nInp, nOut = nOut, nCon = nCon, 
    nShi = nShi, form = form, conDummy = conDummy )
  coefNames <- all.vars( formulaEst )
  coefNames[1] <- "(Intercept)"

  
  #### initialise matrix RMat and vector rVec
  RMat <- matrix( 0, nrow = (nInp + nOut) * nrow( data ), 
    ncol = length( coefNames ) )
  colnames( RMat ) <- coefNames
  rVec <- rep( 0, (nInp + nOut) * nrow( data ) )
  
  # regarding inputs
  for( i in 1:(nInp-1) ) {
    rowsxi <- ( (i-1) * nrow( data ) + 1 ):( i * nrow( data ) )
    for( j in 1:(nInp-1) ) {
      if( i == j ) {
        RMat[ rowsxi, paste0( "alpha_", j ) ] <- 1
      } else{
        RMat[ rowsxi, paste0( "alpha_", j ) ] <- 0
      }
      if( form == "tl" ) {
        for( k in j:(nInp-1) ) {
          if( i == j ) {
            RMat[ rowsxi, paste0( "alpha_", j, "_", k ) ] <-
              log( data[[ xNames[k] ]] / data[[ xNames[nInp] ]] )
          } else if( i == k ) {
            RMat[ rowsxi, paste0( "alpha_", j, "_", k ) ] <-
              log( data[[ xNames[j] ]] / data[[ xNames[nInp] ]] )
          } else{
            RMat[ rowsxi, paste0( "alpha_", j, "_", k ) ] <- 0
          }
        }
        for( k in 1:nOut ) {
          if( i == j ) {
            if( k < nOut ) {
              RMat[ rowsxi, paste0( "psi_", j, "_", k ) ] <-
                datY[[ paste0( "y", k ) ]]
            } else{
              RMat[ rowsxi, paste0( "psi_", j, "_", k ) ] <-
                log( datY$dist_1 )
            }
          } else{
            RMat[ rowsxi, paste0( "psi_", j, "_", k ) ] <- 0
          }
        }
        if( nCon > 0 ) {
          for( k in 1:nCon ) {
            if( i == j ) {
              RMat[ rowsxi, paste0( "xi_", j, "_", k ) ] <-
                data[[ zNames[k] ]]
            } else{
              RMat[ rowsxi, paste0( "xi_", j, "_", k ) ] <- 0
            }
          }
        }
      }
    }
  }
  
  # the "last" input
  rowsxi <- ( (nInp-1) * nrow( data ) + 1 ):( nInp * nrow( data ) )
  for( j in 1:(nInp-1) ) {
    RMat[ rowsxi, paste0( "alpha_", j ) ] <- -1
    if( form == "tl" ) {
      for( k in j:(nInp-1) ) {
        RMat[ rowsxi, paste0( "alpha_", j, "_", k ) ] <-
          -( log( data[[ xNames[j] ]] / data[[ xNames[nInp] ]] ) +
              ( j != k ) * 
              log( data[[ xNames[k] ]] / data[[ xNames[nInp] ]] ) )
      }
      for( k in 1:nOut ) {
        if( k < nOut ) {
          RMat[ rowsxi, paste0( "psi_", j, "_", k ) ] <-
            - datY[[ paste0( "y", k ) ]]
        } else{
          RMat[ rowsxi, paste0( "psi_", j, "_", k ) ] <-
            - log( datY$dist_1 )
        }
      }
      if( nCon > 0 ) {
        for( k in 1:nCon ) {
          RMat[ rowsxi, paste0( "xi_", j, "_", k ) ] <-
            - data[[ zNames[k] ]]
        }
      }
    } 
  }
  rVec[ rowsxi ] <- -1
  
  
  # regarding outputs
  for( i in 1:nOut ) {
    rowsyi <- ( (nInp+i-1)*nrow( data ) + 1 ):( (nInp+i)*nrow( data ) )
    for( j in 1:(nOut-1) ){
      RMat[ rowsyi, paste0( "beta_", j ) ] <-
        -datY[[ paste0( "omega_", j, "_", i ) ]]
    }
    RMat[ rowsyi, paste0( "beta_", nOut ) ] <-
      -data[[ yNames[ i ] ]] / datY$dist_1^2
    if( form == "tl" ) {
      for( j in 1:(nOut-1) ){
        for( k in j:(nOut-1) ){
          RMat[ rowsyi, paste0( "beta_", j, "_", k ) ] <-
            - ( datY[[ paste0( "y", k ) ]] * 
            datY[[ paste0( "omega_", j, "_", i ) ]] +
                ( j!=k ) * datY[[ paste0( "y", j ) ]] * 
                datY[[ paste0( "omega_", k, "_", i ) ]] )
        }
        RMat[ rowsyi, paste0( "beta_", j, "_", nOut ) ] <-
          - datY[[ paste0( "omega_", j, "_", i ) ]] * log( datY$dist_1 ) -
          datY[[ paste0( "y", j ) ]] * data[[ yNames[ i ] ]] / datY$dist_1^2
      }
      RMat[ rowsyi, paste0( "beta_", nOut, "_", nOut ) ] <-
        - log( datY$dist_1 ) * data[[ yNames[ i ] ]] / datY$dist_1^2
      for( j in 1:(nInp-1) ){
        for( k in 1:(nOut-1) ){
          RMat[ rowsyi, paste0( "psi_", j, "_", k ) ] <-
           - log( data[[ xNames[j] ]] / data[[ xNames[nInp] ]] ) * 
            datY[[ paste0( "omega_", k, "_", i ) ]]
        }
        RMat[ rowsyi, paste0( "psi_", j, "_", nOut ) ] <-
          - log( data[[ xNames[j] ]] / data[[ xNames[nInp] ]] ) * 
          data[[ yNames[ i ] ]] / datY$dist_1^2
      }
      if( nCon > 0 ) {
        for( k in 1:(nCon) ){
          for( j in 1:(nOut-1) ){
            RMat[ rowsyi, paste0( "zeta_", j, "_", k ) ] <-
              -data[[ zNames[k] ]] * datY[[ paste0( "omega_", j, "_", i ) ]]
          }
          RMat[ rowsyi, paste0( "zeta_", nOut, "_", k ) ] <-
            -data[[ zNames[k] ]] * data[[ yNames[ i ] ]] / datY$dist_1^2
        }
      }
    }
  }
  
  result <- list( RMat = RMat, rVec = rVec )
  return( result )
}