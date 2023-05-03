library( "micEconDistRay" )
library( "quadprog" )

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
## IGNORE_RDIFF_BEGIN
cbind( names( estCD ) )
## IGNORE_RDIFF_END
lapply( estCD$coefList, function(x) round( x, 3 ) )
apply( estCD$mono, 2, table )
lapply( estCD$ela, function(x) round( summary(x), 2 ) )

# calculate the dependent variable (logarithm of predicted distance)
appleProdFr86$distCD <- distRayCalc( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86, 
  coef = coef( estCD ), form = "cd" )
round( summary( appleProdFr86$distCD ), 3 )

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

# vector of unrestricted coefficients and their covariance matrix
nCoefCD <- length( coef( estCD ) ) - 2
uCoefCD <- coef( estCD )[ 1:nCoefCD ]
uCovInvCD <- solve( vcov( estCD )[ 1:nCoefCD, 1:nCoefCD ] )

# obtain the matrix and vector to impose monotonicity
restrCD <- distRayMonoRestr( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86,
  form = "cd" )

# obtain the restricted coefficients
minDistCD <- solve.QP( Dmat = uCovInvCD, dvec = rep( 0, nCoefCD ),
  Amat = t( restrCD$RMat ), bvec = - restrCD$RMat %*% uCoefCD + restrCD$rVec )
rCoefCD <- minDistCD$solution + uCoefCD
round( rCoefCD, 3 )

# calculate elasticities based on restricted coefficients
rElaCD <- distRayEla( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = rCoefCD, data = appleProdFr86, form = "cd" )
lapply( rElaCD, function(x) round( summary(x), 2 ) )


### Translog ray-based input distance function
estTL <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86 )
cbind( round( coef( estTL ), 2 ) )
## IGNORE_RDIFF_BEGIN
cbind( names( estTL ) )
## IGNORE_RDIFF_END
lapply( estTL$coefList, function(x) round( x, 3 ) )
apply( estTL$mono, 2, table )
lapply( estTL$ela, function(x) round( summary(x), 2 ) )

# calculate the dependent variable (logarithm of predicted distance)
appleProdFr86$logDistTL <- distRayCalc( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86, 
  coef = coef( estTL ) )
round( summary( appleProdFr86$logDistTL ), 3 )

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

# vector of unrestricted coefficients and their covariance matrix
nCoefTL <- length( coef( estTL ) ) - 2
uCoefTL <- coef( estTL )[ 1:nCoefTL ]
uCovInvTL <- solve( vcov( estTL )[ 1:nCoefTL, 1:nCoefTL ] )

# obtain the matrix and vector to impose monotonicity
restrTL <- distRayMonoRestr( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86 )

# obtain the restricted coefficientslibrary( "quadprog" )
minDistTL <- solve.QP( Dmat = uCovInvTL, dvec = rep( 0, nCoefTL ),
  Amat = t( restrTL$RMat ), bvec = - restrTL$RMat %*% uCoefTL + restrTL$rVec )
rCoefTL <- minDistTL$solution + uCoefTL
round( rCoefTL, 3 )

# calculate elasticities based on restricted coefficients
rElaTL <- distRayEla( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = rCoefTL, data = appleProdFr86 )
lapply( rElaTL, function(x) round( summary(x), 2 ) )

