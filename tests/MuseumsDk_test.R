# a part of this test script is a partial replication of Price & Henningsen
# (forthcoming): A Ray-Based Input Distance Function to Model Zero-Valued
# Output Quantities: Derivation and an Empirical Application, Journal of
# Productivity Analysis.

library( "micEconDistRay" )
library( "quadprog" )

data( "MuseumsDk" )

# prepare variables
MuseumsDk$pub <- with( MuseumsDk, aarc + ach + aah + anh ) 
MuseumsDk$k <- with( MuseumsDk, expProperty / ipc )
MuseumsDk$fte <- MuseumsDk$ftesc + MuseumsDk$ftensc
MuseumsDk$logUnits <- log( MuseumsDk$units )

# inputs
xNames <- c( "k", "ftesc", "ftensc" )
# outputs 
yNames <- c( "vis", "expCons", "pub", "exh", "edu", "ev" )
# control variables 
zNames <- c( "logUnits", "resp" )
conDummy <- c( 2 )

# remove observations that cannot be used in the estimations
MuseumsDk <- subset( MuseumsDk, rowSums( MuseumsDk[ , xNames ] <= 0 ) == 0 &
    rowSums( is.na( MuseumsDk[ , xNames ] ) ) == 0 &
    rowSums( is.na( MuseumsDk[ , yNames ] ) ) == 0 )

# remove observation with no visitors
MuseumsDk <- subset( MuseumsDk, vis > 0 )

# mean-scaling input and output quantities
for( n in c( xNames, yNames ) ) {
  MuseumsDk[[ n ]] <- MuseumsDk[[ n ]] / mean( MuseumsDk[[ n ]] ) 
}


### estimate a ray-based Cobb-Douglas input distance function
estCD <- distRayEst( xNames = xNames, yNames = yNames, zNames = zNames,
  data = MuseumsDk, form = "cd" )
cbind( round( coef( estCD ), 2 ) )
cbind( names( estCD ) )
lapply( estCD$coefList, function(x) round( x, 3 ) )
apply( estCD$mono, 2, table )
lapply( estCD$ela, function(x) round( summary(x), 2 ) )

# calculate the dependent variable (logarithm of predicted distance)
MuseumsDk$distCD <- distRayCalc( 
  xNames = xNames, yNames = yNames, zNames = zNames, 
  data = MuseumsDk, coef = coef( estCD ), form = "cd" )
round( summary( MuseumsDk$distCD ), 3 )

# calculate elasticities
elaCD <- distRayEla( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = coef( estCD ), data = MuseumsDk, form = "cd" )
all.equal( elaCD, estCD$ela )
lapply( elaCD, function(x) round( summary(x), 2 ) )

# calculate derivatives
derivCD <- distRayDeriv( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = coef( estCD ), data = MuseumsDk, form = "cd" )
lapply( derivCD, function(x) round( summary(x), 
  3 - round( log( max( abs( x ) ), 10 ) ) ) )

# vector of unrestricted coefficients and their covariance matrix
nCoefCD <- length( coef( estCD ) ) - 2
uCoefCD <- coef( estCD )[ 1:nCoefCD ]
uCovInvCD <- solve( vcov( estCD )[ 1:nCoefCD, 1:nCoefCD ] )

# obtain the matrix and vector to impose monotonicity
restrCD <- distRayMonoRestr( 
  xNames = xNames, yNames = yNames, zNames = zNames, 
  data = MuseumsDk, form = "cd" )

# obtain the restricted coefficients
minDistCD <- solve.QP( Dmat = uCovInvCD, dvec = rep( 0, nCoefCD ),
  Amat = t( restrCD$RMat ), bvec = - restrCD$RMat %*% uCoefCD + restrCD$rVec )
rCoefCD <- minDistCD$solution + uCoefCD
round( rCoefCD, 3 )

# calculate elasticities based on restricted coefficients
rElaCD <- distRayEla( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = rCoefCD, data = MuseumsDk, form = "cd" )
lapply( rElaCD, function(x) round( summary(x), 2 ) )


### estimate a ray-based Cobb-Douglas input distance function
### with control variables as "shifters" only
estTLShift <- distRayEst( xNames = xNames, yNames = yNames, sNames = zNames,
  data = MuseumsDk, form = "tl" )
cbind( round( coef( estTLShift ), 2 ) )
cbind( names( estTLShift ) )
lapply( estTLShift$coefList, function(x) round( x, 3 ) )
apply( estTLShift$mono, 2, table )
lapply( estTLShift$ela, function(x) round( summary(x), 2 ) )

# calculate the dependent variable (logarithm of predicted distance)
MuseumsDk$distTLShift <- distRayCalc( 
  xNames = xNames, yNames = yNames, sNames = zNames, 
  data = MuseumsDk, coef = coef( estTLShift ) )
round( summary( MuseumsDk$distTLShift ), 3 )

# calculate elasticities
elaTLShift <- distRayEla( xNames = xNames, yNames = yNames, sNames = zNames,
  coef = coef( estTLShift ), data = MuseumsDk )
all.equal( elaTLShift, estTLShift$ela )
lapply( elaTLShift, function(x) round( summary(x), 2 ) )

# calculate derivatives
derivTLShift <- distRayDeriv( xNames = xNames, yNames = yNames, sNames = zNames,
  coef = coef( estTLShift ), data = MuseumsDk )
lapply( derivTLShift, function(x) round( summary(x), 
  3 - round( log( max( abs( x ) ), 10 ) ) ) )

# vector of unrestricted coefficients and their covariance matrix
nCoefTLShift <- length( coef( estTLShift ) ) - 2
uCoefTLShift <- coef( estTLShift )[ 1:nCoefTLShift ]
uCovInvTLShift <- solve( vcov( estTLShift )[ 1:nCoefTLShift, 1:nCoefTLShift ] )

# obtain the matrix and vector to impose monotonicity
restrTLShift <- distRayMonoRestr( 
  xNames = xNames, yNames = yNames, sNames = zNames, 
  data = MuseumsDk )

# obtain the restricted coefficients
minDistTLShift <- solve.QP( Dmat = uCovInvTLShift, dvec = rep( 0, nCoefTLShift ),
  Amat = t( restrTLShift$RMat ), bvec = - restrTLShift$RMat %*% uCoefTLShift + restrTLShift$rVec )
rCoefTLShift <- minDistTLShift$solution + uCoefTLShift
round( rCoefTLShift, 3 )

# calculate elasticities based on restricted coefficients
rElaTLShift <- distRayEla( xNames = xNames, yNames = yNames, sNames = zNames,
  coef = rCoefTLShift, data = MuseumsDk )
lapply( rElaTLShift, function(x) round( summary(x), 2 ) )


### estimate a ray-based "full" Translog input distance function
estTLFull <- distRayEst( xNames = xNames, yNames = yNames, zNames = zNames,
  data = MuseumsDk, form = "tl" )
cbind( round( coef( estTLFull ), 2 ) )
cbind( names( estTLFull ) )
lapply( estTLFull$coefList, function(x) round( x, 3 ) )
apply( estTLFull$mono, 2, table )
lapply( estTLFull$ela, function(x) round( summary(x), 2 ) )

# calculate the dependent variable (logarithm of predicted distance)
MuseumsDk$distTLFull <- distRayCalc( 
  xNames = xNames, yNames = yNames, zNames = zNames, 
  data = MuseumsDk, coef = coef( estTLFull ) )
round( summary( MuseumsDk$distTLFull ), 3 )

# calculate elasticities
elaTLFull <- distRayEla( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = coef( estTLFull ), data = MuseumsDk, conDummy = conDummy )
all.equal( elaTLFull, estTLFull$ela )
lapply( elaTLFull, function(x) round( summary(x), 2 ) )

# calculate derivatives
derivTLFull <- distRayDeriv( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = coef( estTLFull ), data = MuseumsDk )
lapply( derivTLFull, function(x) round( summary(x), 
  3 - round( log( max( abs( x ) ), 10 ) ) ) )

# vector of unrestricted coefficients and their covariance matrix
nCoefTLFull <- length( coef( estTLFull ) ) - 2
uCoefTLFull <- coef( estTLFull )[ 1:nCoefTLFull ]
uCovInvTLFull <- solve( vcov( estTLFull )[ 1:nCoefTLFull, 1:nCoefTLFull ] )

# obtain the matrix and vector to impose monotonicity
restrTLFull <- distRayMonoRestr( 
  xNames = xNames, yNames = yNames, zNames = zNames, 
  data = MuseumsDk, conDummy = conDummy )

# obtain the restricted coefficients
minDistTLFull <- solve.QP( Dmat = uCovInvTLFull, dvec = rep( 0, nCoefTLFull ),
  Amat = t( restrTLFull$RMat ), bvec = - restrTLFull$RMat %*% uCoefTLFull + restrTLFull$rVec )
rCoefTLFull <- minDistTLFull$solution + uCoefTLFull
round( rCoefTLFull, 3 )

# calculate elasticities based on restricted coefficients
rElaTLFull <- distRayEla( xNames = xNames, yNames = yNames, zNames = zNames,
  coef = rCoefTLFull, data = MuseumsDk, conDummy = conDummy )
lapply( rElaTLFull, function(x) round( summary(x), 2 ) )

