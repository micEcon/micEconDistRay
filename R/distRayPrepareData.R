distRayPrepareData <- function( xNames, yNames, 
  zNames = NULL, sNames = NULL, eNames = NULL,
  data, form, fixThetas = FALSE ) {

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
  nEff <- length( eNames )
  
  
  #### calculate angles and distances of output quantities ####
  datY <- distRayCalcYVar( yNames, data, calcOmegas = TRUE, 
    fixThetas = fixThetas )
  
  
  #### initialise data set that should be returned ####
  datEst <- data.frame( row.names = row.names( data ) )
  
  # dependent variable
  datEst$dep <- - log( data[[ xNames[nInp] ]] )
  
  # alpha_i
  for( i in 1:( nInp - 1 ) ) {
    datEst[[ paste0( "alpha_", i ) ]] <- 
      log( data[[ xNames[i] ]] / data[[ xNames[nInp] ]] )
  }

  # beta_i 
  for( i in 1:nOut ) {
    datEst[[ paste0( "beta_", i ) ]] <- datY[[ paste0( "y", i ) ]]
  }

  # delta_i
  if( nCon > 0 ) {
    for( i in 1:nCon ) {
      datEst[[ paste0( "delta_", i ) ]] <- data[[ zNames[i] ]]
    }
  }
  if( nShi > 0 ) {
    for( i in 1:nShi ) {
      datEst[[ paste0( "delta_", nCon + i ) ]] <- data[[ sNames[i] ]]
    }
  }
  
  if( form == "tl" ) {
    # alpha_ij
    for( i in 1:( nInp - 1 ) ) {
      for( j in i:( nInp - 1 ) ) {
        datEst[[ paste0( "alpha_", i, "_", j ) ]] <- ifelse( i == j, 0.5, 1 ) * 
          log( data[[ xNames[i] ]] / data[[ xNames[nInp] ]] ) *
          log( data[[ xNames[j] ]] / data[[ xNames[nInp] ]] )
      }  
    }
    
    # beta_ij 
    for( i in 1:nOut ) {
      for( j in i:nOut ) {
        datEst[[ paste0( "beta_", i, "_", j ) ]] <- ifelse( i == j, 0.5, 1 ) *
          datY[[ paste0( "y", i ) ]] * datY[[ paste0( "y", j ) ]]
      }
    }
    
    # psi_ij
    for( i in 1:( nInp - 1 ) ) {
      for( j in 1:nOut ) {
        datEst[[ paste0( "psi_", i, "_", j ) ]] <- 
          log( data[[ xNames[i] ]] / data[[ xNames[nInp] ]] ) *
          datY[[ paste0( "y", j ) ]] 
      }
    }
    
    if( nCon > 0 ) {
      # delta_ij
      for( i in 1:nCon ) {
        for( j in i:nCon ) {
          datEst[[ paste0( "delta_", i, "_", j ) ]] <- ifelse( i == j, 0.5, 1 ) *
            data[[ zNames[i] ]]  * data[[ zNames[j] ]]
        }
      }
      
      # xi_ij
      for( i in 1:( nInp - 1 ) ) {
        for( j in 1:nCon ) {
          datEst[[ paste0( "xi_", i, "_", j ) ]] <-
            log( data[[ xNames[i] ]] / data[[ xNames[nInp] ]] ) *
            data[[ zNames[j] ]]
        }
      }

      # zeta_ij    
      for( i in 1:nOut ) {
        for( j in 1:nCon ) {
          datEst[[ paste0( "zeta_", i, "_", j ) ]] <-
            datY[[ paste0( "y", i ) ]] * data[[ zNames[j] ]]
        }
      }
    }
  }
  
  if( nEff > 0 ) {
    # phi_i 
    for( i in 1:nEff ) {
      datEst[[ paste0( "phi_", i ) ]] <- data[[ eNames[i] ]]
    }
  }

  return( datEst )
}