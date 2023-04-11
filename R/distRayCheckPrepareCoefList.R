distRayCheckPrepareCoefList <- function( nInp, nOut, nCon, nShi, coef, form,
  conDummy = NULL ) {

  # create 'coefList'
  if( is.list( coef ) ) {
    coefList <- coef
  } else if( is.vector( coef ) && is.numeric( coef ) ) {
    coefList <- distRayPrepareCoef( nInp = nInp, nOut = nOut, nCon = nCon, 
      nShi = nShi, coefVec = coef, form = form, conDummy = NULL )
    if( form == "tl" && nCon > 0 ) {
      for( i in 1:nCon ) {
        if( is.na( coefList$deltaMat[ i, i ] ) && i %in% conDummy ) {
          coefList$deltaMat[ i, i ] <- 0
        }
      }
    }
  } else {
    return( "argument 'coef' must be a list or a numeric vector" )
  }
  
  #### check argument coefList ####
  #### and give an informative error message if something is incorrect
  if( !inherits( coefList, "list" ) ) {
    return( "argument 'coefList' must be a list" )
  }
  expectedCoef <- c( "alpha0", "alphaVec", "betaVec" )
  if( ( nCon + nShi ) > 0 ) {
    expectedCoef <- c( expectedCoef, "deltaVec" )
  }
  if( form == "tl" ) { 
    expectedCoef <- c( expectedCoef, "alphaMat", "betaMat", "psiMat" )
    if( nCon > 0 ) {
      expectedCoef <- c( expectedCoef, "deltaMat", "xiMat", "zetaMat" )
    }
  }
  unexpectedCoef <- names( coefList )[ !names( coefList ) %in% expectedCoef ] 
  if( length( unexpectedCoef ) > 0 ) {
    return( paste( "argument 'coefList' has unexpected components:",
      paste( unexpectedCoef, collapse = ", " ) ) )
  }
  missingCoef <- expectedCoef[ !expectedCoef %in% names( coefList ) ] 
  if( length( missingCoef ) > 0 ) {
    return( paste( "argument 'coefList' has missing components:",
      paste( missingCoef, collapse = ", " ) ) )
  }
  if( length( coefList$alpha0 ) != 1 ||
      !is.numeric( coefList$alpha0 ) ) {
    return( paste( "component 'alpha0' of argument 'coefList' must be a",
      "numeric scalar" ) )
  }
  if( length( coefList$alphaVec ) != nInp || 
      !is.vector( coefList$alphaVec ) || 
      !is.numeric( coefList$alphaVec ) ) {
    return( paste( "component 'alphaVec' of argument 'coefList' must be a",
      "numeric vector of length", nInp ) )
  }
  if( abs( sum( coefList$alphaVec ) - 1 ) > 1e-6 ) {
    return( paste( "the elements of component 'alphaVec'",
      "of argument 'coefList' must sum up to one" ) )
  }
  if( length( coefList$betaVec ) != nOut || 
      !is.vector( coefList$betaVec ) || 
      !is.numeric( coefList$betaVec ) ) {
    return( paste( "component 'betaVec' of argument 'coefList' must be a",
      "numeric vector of length", nOut ) )
  }
  if( ( nCon + nShi ) > 0 ) {
    if( length( coefList$deltaVec ) != ( nCon + nShi ) || 
        !is.vector( coefList$deltaVec ) || 
        !is.numeric( coefList$deltaVec ) ) {
      return( paste( "component 'deltaVec' of argument 'coefList' must be a",
        "numeric vector of length", nCon + nShi ) )
    }
  }
  if( form == "tl" ) {
    if( nrow( coefList$alphaMat ) != nInp ||
        ncol( coefList$alphaMat ) != nInp ||
        !is.matrix( coefList$alphaMat ) || 
        !is.numeric( coefList$alphaMat ) ) {
      return( paste( "component 'alphaMat' of argument 'coefList' must be a",
        "numeric matrix of dimension", nInp, "x", nInp ) )
    }
    if( !isSymmetric( coefList$alphaMat ) ) {
      return( paste( "component 'alphaMat' of argument 'coefList'",
        "must be a symmetric matrix" ) )
    }
    if( max( abs( c( rowSums( coefList$alphaMat ), 
      colSums( coefList$alphaMat ) ) ) ) > 1e-6 ) {
      return( paste( "all rows and columns of component 'alphaMat'",
        "of argument 'coefList' must sum up to zero" ) )
    }
    if( nrow( coefList$betaMat ) != nOut ||
        ncol( coefList$betaMat ) != nOut ||
        !is.matrix( coefList$betaMat ) || 
        !is.numeric( coefList$betaMat ) ) {
      return( paste( "component 'betaMat' of argument 'coefList' must be a",
        "numeric matrix of dimension", nOut, "x", nOut ) )
    }
    if( !isSymmetric( coefList$betaMat ) ) {
      return( paste( "component 'betaMat' of argument 'coefList'",
        "must be a symmetric matrix" ) )
    }
    if( nrow( coefList$psiMat ) != nInp ||
        ncol( coefList$psiMat ) != nOut ||
        !is.matrix( coefList$psiMat ) || 
        !is.numeric( coefList$psiMat ) ) {
      return( paste( "component 'psiMat' of argument 'coefList' must be a",
        "numeric matrix of dimension", nInp, "x", nOut ) )
    }
    if( max( abs( colSums( coefList$psiMat ) ) ) > 1e-6 ) {
      return( paste( "all columns of component 'psiMat'",
        "of argument 'coefList' must sum up to zero" ) )
    }
    if( nCon > 0 ) {
      if( nrow( coefList$deltaMat ) != nCon ||
          ncol( coefList$deltaMat ) != nCon ||
          !is.matrix( coefList$deltaMat ) || 
          !is.numeric( coefList$deltaMat ) ) {
        return( paste( "component 'deltaMat' of argument 'coefList' must be a",
          "numeric matrix of dimension", nCon, "x", nCon ) )
      }
      if( !isSymmetric( coefList$deltaMat ) ) {
        return( paste( "component 'deltaMat' of argument 'coefList'",
          "must be a symmetric matrix" ) )
      }
      if( nrow( coefList$xiMat ) != nInp ||
          ncol( coefList$xiMat ) != nCon ||
          !is.matrix( coefList$xiMat ) || 
          !is.numeric( coefList$xiMat ) ) {
        return( paste( "component 'xiMat' of argument 'coefList' must be a",
          "numeric matrix of dimension", nInp, "x", nCon ) )
      }
      if( max( abs( colSums( coefList$xiMat ) ) ) > 1e-6 ) {
        return( paste( "all columns of component 'xiMat'",
          "of argument 'coefList' must sum up to zero" ) )
      }
      if( nrow( coefList$zetaMat ) != nOut ||
          ncol( coefList$zetaMat ) != nCon ||
          !is.matrix( coefList$zetaMat ) || 
          !is.numeric( coefList$zetaMat ) ) {
        return( paste( "component 'zetaMat' of argument 'coefList' must be a",
          "numeric matrix of dimension", nOut, "x", nCon ) )
      }
    }
  }
  return( coefList )
}