library( "micEconDistRay" )

# load and prepare data set 
data( germanFarms, package = "micEcon" )
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

# Cobb-Douglas ray-based input distance function
estCD <- distRayEst( xNames = c( "land", "qLabor", "qVarInput" ),
  yNames = c( "vCrop", "vAnimal" ),
  data = germanFarms, form = "cd" )
cbind( round( coef( estCD ), 2 ) )
cbind( names( estCD ) )
lapply( estCD$coefList, function(x) round( x, 3 ) )
apply( estCD$mono, 2, table )
lapply( estCD$ela, function(x) round( summary(x), 2 ) )

# Translog ray-based input distance function
estTL <- suppressWarnings( distRayEst( xNames = c( "land", "qLabor", "qVarInput" ),
  yNames = c( "vCrop", "vAnimal" ),
  data = germanFarms ) )
cbind( round( coef( estTL ), 2 ) )
cbind( names( estTL ) )
lapply( estTL$coefList, function(x) round( x, 3 ) )
apply( estTL$mono, 2, table )
lapply( estTL$ela, function(x) round( summary(x), 2 ) )
