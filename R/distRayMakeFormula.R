distRayMakeFormula <- function( nInp, nOut, nCon = 0, nShi = 0, nEff = 0,
  form, conDummy = NULL ) {

  #### check arguments nInp, nOut, nCon, nShi, conDummy ####
  testRes <- distRayCheckNumVars( nInp = nInp, nOut = nOut, nCon = nCon, 
    nShi = nShi, form = form, conDummy = conDummy )
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
  
  #### create a character string that defines the regression formula ####
  estForm <- "dep ~ 1"

  # alpha_i
  if( nInp > 1 ) {
    estForm <- paste( estForm, "+",
      paste0( "alpha_", 1:( nInp - 1 ), collapse = " + " ) )
  }

  # beta_i
  estForm <- paste( estForm, "+", paste0( "beta_", 1:nOut, collapse = " + " ) )

  # delta_i
  if( nCon + nShi > 0 ) {
    estForm <- paste( estForm, "+", paste0( "delta_", 1:( nCon + nShi ), 
      collapse = " + " ) )
  }
  
  if( form == "tl" ) {
    # alpha_ij
    if( nInp > 1 ) {
      for( i in 1:(nInp-1) ) {
        for( j in i:(nInp-1) ) {
          estForm <- paste0( estForm, " + alpha_", i, "_", j )
        }
      }
    }
    
    # beta_ij
    for( i in 1:nOut ) {
      for( j in i:nOut ) {
        estForm <- paste0( estForm, " + beta_", i, "_", j )
      } 
    }
    
    # psi_ij
    for( i in 1:( nInp - 1 ) ) {
      for( j in 1:nOut ) {
        estForm <- paste0( estForm, " + psi_", i, "_", j )
      }
    }
    
    if( nCon > 0 ) {
      # delta_ij
      for( i in 1:nCon ) {
        for( j in i:nCon ) {
          if( i != j || ! i %in% conDummy ) {
            estForm <- paste0( estForm, " + delta_", i, "_", j )
          }
        }
      }

      # xi_ij
      if( nInp > 1 ) {
        for( i in 1:( nInp - 1 ) ) {
          for( j in 1:nCon ) {
            estForm <- paste0( estForm, " + xi_", i, "_", j )
          }
        }
      }
    
      # zeta_ij
      for( i in 1:nOut ) {
        for( j in 1:nCon ) {
          estForm <- paste0( estForm, " + zeta_", i, "_", j )
        }
      }
    }
  }
  
  # phi_i
  if( nEff > 0 ) {
    estForm <- paste( estForm, "|", paste0( "phi_", 1:nEff, 
      collapse = " + " ), "- 1" )
  }
  
  # convert the character string to a formula object
  estForm <- as.formula( estForm )
  
  return( estForm )
}
