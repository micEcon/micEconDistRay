library( "micEconDistRay" )

# load and prepare data set 
data( appleProdFr86, package = "micEcon" )
appleProdFr86$qCap <- appleProdFr86$vCap / appleProdFr86$pCap
appleProdFr86$qLab <- appleProdFr86$vLab / appleProdFr86$pLab
appleProdFr86$qMat <- appleProdFr86$vMat / appleProdFr86$pMat

### Cobb-Douglas ray-based input distance function
estCD <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86, form = "cd" )
cbind( round( coef( estCD ), 2 ) )
cbind( names( estCD ) )
lapply( estCD$coefList, function(x) round( x, 3 ) )
apply( estCD$mono, 2, table )
lapply( estCD$ela, function(x) round( summary(x), 2 ) )

# calculate elasticities
elaCD <- distRayEla( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estCD ), data = appleProdFr86, form = "cd" )
all.equal( elaCD, estCD$ela )
lapply( elaCD, function(x) round( summary(x), 2 ) )

# calculate derivatives
derivCD <- distRayDeriv( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estCD ), data = appleProdFr86, form = "cd" )
lapply( derivCD, function(x) round( summary(x), 
  3 - round( log( max( abs( x ) ), 10 ) ) ) )


### Translog ray-based input distance function
estTL <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86 )
cbind( round( coef( estTL ), 2 ) )
cbind( names( estTL ) )
lapply( estTL$coefList, function(x) round( x, 3 ) )
apply( estTL$mono, 2, table )
lapply( estTL$ela, function(x) round( summary(x), 2 ) )

# calculate elasticities
elaTL <- distRayEla( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estTL ), data = appleProdFr86 )
all.equal( elaTL, estTL$ela )
lapply( elaTL, function(x) round( summary(x), 2 ) )

# calculate derivatives
derivTL <- distRayDeriv( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estTL ), data = appleProdFr86 )
lapply( derivTL, function(x) round( summary(x), 
  3 - round( log( max( abs( x ) ), 10 ) ) ) )

