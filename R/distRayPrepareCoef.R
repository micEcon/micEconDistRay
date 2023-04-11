distRayPrepareCoef <- function( nInp, nOut, nCon = 0, nShi = 0, coefVec, 
  form, conDummy = NULL ) {

  #### check arguments nInp, nOut, nCon, nShi, conDummy ####
  testRes <- distRayCheckNumVars( nInp = nInp, nOut = nOut, nCon = nCon, 
    nShi = nShi, form = form, conDummy = conDummy )
  if( !is.null( testRes ) ) {
    stop( testRes )
  }

  #### check argument 'coefVec' ####
  if( !is.vector( coefVec) || !is.numeric( coefVec ) ) {
    stop( "argument 'coefVec' must be a vector of numeric values" )
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
  
  # create list for the different components of the coefficients 
  coefList <- list()
  
  # alpha_0
  coefList$alpha0 <- coefVec[ "(Intercept)"]
  
  # alpha_i coefficients
  if( nInp > 1 ) {
    coefList$alphaVec <- unname( coefVec[ paste0( "alpha_", 1:( nInp - 1 ) ) ] )
    coefList$alphaVec <- c( coefList$alphaVec, 1 - sum( coefList$alphaVec ) )
  }
  
  # beta_i coefficients
  coefList$betaVec <- unname( coefVec[ paste0( "beta_", 1:nOut ) ] )

  # delta_i coefficients
  if( nCon + nShi > 1 ) {
    coefList$deltaVec <- unname( coefVec[ paste0( "delta_", 1:( nCon + nShi ) ) ] )
  }

  if( form == "tl" ) {
    # alpha_i_j coefficients
    if( nInp > 1 ) {
      coefList$alphaMat <- matrix( NA, nrow = nInp - 1, ncol = nInp - 1 )
      for( i in 1:(nInp-1) ){
        for( j in 1:(nInp-1) ) {
          coefList$alphaMat[i,j] <- coefVec[ 
            paste0( "alpha_", min( i, j ), "_", max( i, j ) ) ]
        }
      }
      coefList$alphaMat <- rbind( coefList$alphaMat, -colSums( coefList$alphaMat ) )
      coefList$alphaMat <- cbind( coefList$alphaMat, -rowSums( coefList$alphaMat ) )
    }
    
    # beta_i_j coefficients 
    coefList$betaMat <- matrix( NA, nrow = nOut, ncol = nOut )
    for( i in 1:(nOut) ){
      for( j in 1:(nOut) ) {
        coefList$betaMat[i,j] <- coefVec[ 
          paste0( "beta_", min( i, j ), "_", max( i, j ) ) ]
      }
    }
    
    # delta_i_j coefficients
    if( nCon > 0 ) {
      coefList$deltaMat <- matrix( NA, nrow = nCon, ncol = nCon )
      for( i in 1:(nCon) ){
        for( j in 1:(nCon) ) {
          coefList$deltaMat[i,j] <- coefVec[ 
            paste0( "delta_", min( i, j ), "_", max( i, j ) ) ]
        }
      }
      if( !is.null( conDummy) ) {
        for( i in conDummy ){
          if( is.na( coefList$deltaMat[i,i] ) ) {
            coefList$deltaMat[i,i] <- 0
          } else {
            stop( "control variable ", i, " cannot be a dummy variable,",
               " because coefficient delta_{", i, ",", i, "} has bee estimated" )
          }
        }
      }
    }
    
    # psi_i_j coefficients
    if( nInp > 1 ) {
      coefList$psiMat <- matrix( NA, nrow = nInp - 1, ncol = nOut )
      for( i in 1:(nInp-1) ){
        for( j in 1:(nOut) ) {
          coefList$psiMat[i,j] <- coefVec[ 
            paste0( "psi_", i, "_", j ) ]
        }
      }
      coefList$psiMat <- rbind( coefList$psiMat, -colSums( coefList$psiMat ) )
    }
    
    # xi_i_j coefficients
    if( nInp > 1 && nCon > 0 ) {
      coefList$xiMat <- matrix( NA, nrow = nInp - 1, ncol = nCon )
      for( i in 1:(nInp-1) ){
        for( j in 1:(nCon) ) {
          coefList$xiMat[i,j] <- coefVec[ 
            paste0( "xi_", i, "_", j ) ]
        }
      }
      coefList$xiMat <- rbind( coefList$xiMat, -colSums( coefList$xiMat ) )
    }
    
    # zeta_i_j coefficients
    if( nCon > 0 ) {
      coefList$zetaMat <- matrix( NA, nrow = nOut, ncol = nCon )
      for( i in 1:(nOut) ){
        for( j in 1:(nCon) ) {
          coefList$zetaMat[i,j] <- coefVec[ 
            paste0( "zeta_", i, "_", j ) ]
        }
      }
    }
  }
  
  return( coefList )
}
