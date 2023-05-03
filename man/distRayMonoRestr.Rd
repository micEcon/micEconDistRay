\name{distRayMonoRestr}
\alias{distRayMonoRestr}
\title{Imposing Monotonicity on a Ray-Based Input Distance Function}

\description{
   Create a matrix and vector for imposing monotonicity 
   on a ray-based input distance function
   (Price & Henningsen, forthcomong).
}

\usage{
distRayMonoRestr( xNames, yNames, zNames = NULL, sNames = NULL, 
  data, form = "tl", conDummy = NULL )
}

\arguments{
  \item{xNames}{a vector of character strings 
    containing the names of the variables
    that indicate the input quantities.}
  \item{yNames}{a vector of two or more character strings 
    containing the names of the variables 
    that indicate the output quantities.}
  \item{zNames}{an optional vector of character strings 
    containing the names of \sQuote{environmental} variables,
    i.e., variables that affect the production possibility set
    (i.e., the feasible combinations of input-output quantities)
    that---in the case of a Translog functional form---should be interacted 
    with the input quantities and the angles of the output vector.}
  \item{sNames}{an optional vector of character strings 
    containing the names of \sQuote{environmental} variables,
    i.e., variables that affect the production possibility set
    (i.e., the feasible combinations of input-output quantities)
    that---in the case of a Translog functional form---should not be interacted 
    with the input quantities and the angles of the output vector.}
  \item{data}{data frame containing the data.}
  \item{form}{a character string that indicates the functional form;
    currently, \code{"cd"} for the Cobb-Douglas functional form
    and \code{"tl"} for the Translog functional form are available.}
  \item{conDummy}{an optional numeric vector indicating the positions
    in argument \code{zNames} that are dummy variables
    so that quadratic terms of these variables are omitted.}
}

\value{
   A list that contains a matrix (\code{RMat}) 
   and a vector (\code{rVec}).
}

\author{Arne Henningsen and Juan José Price}

\references{
  Price, J.J. & Henningsen, A. (forthcoming):
  A Ray-Based Input Distance Function to Model Zero-Valued Output Quantities: 
  Derivation and an Empirical Application.
  Journal of Productivity Analysis.
}

\examples{
# load and prepare data set 
data( appleProdFr86, package = "micEcon" )
appleProdFr86$qCap <- appleProdFr86$vCap / appleProdFr86$pCap
appleProdFr86$qLab <- appleProdFr86$vLab / appleProdFr86$pLab
appleProdFr86$qMat <- appleProdFr86$vMat / appleProdFr86$pMat

# Cobb-Douglas ray-based input distance function (with manually set parameters)
estCD <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86, form = "cd" )
summary( estCD )

# the vector of unrestricted coefficients and their covariance matrix
nCoefCD <- length( coef( estCD ) ) - 2
uCoefCD <- coef( estCD )[ 1:nCoefCD ]
uCovInvCD <- solve( vcov( estCD )[ 1:nCoefCD, 1:nCoefCD ] )

# obtain the matrix and vector to impose monotonicity
restrCD <- distRayMonoRestr( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86,
  form = "cd" )

# obtain the restricted coefficients
library( "quadprog" )
minDistCD <- solve.QP( Dmat = uCovInvCD, dvec = rep( 0, nCoefCD ),
  Amat = t( restrCD$RMat ), bvec = - restrCD$RMat \%*\% uCoefCD + restrCD$rVec )
rCoefCD <- minDistCD$solution + uCoefCD


# Translog ray-based input distance function (with estimated parameters)
estTL <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86 )
appleProdFr86$logDistTL <- distRayCalc( xNames = c( "qCap", "qLab", "qMat" ),
  yNames= c( "qApples", "qOtherOut" ), data = appleProdFr86,
  coef = coef( estTL ) )
summary( appleProdFr86$logDistTL )

# the vector of unrestricted coefficients and their covariance matrix
nCoefTL <- length( coef( estTL ) ) - 2
uCoefTL <- coef( estTL )[ 1:nCoefTL ]
uCovInvTL <- solve( vcov( estTL )[ 1:nCoefTL, 1:nCoefTL ] )

# obtain the matrix and vector to impose monotonicity
restrTL <- distRayMonoRestr( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86 )

# obtain the restricted coefficientslibrary( "quadprog" )
minDistTL <- solve.QP( Dmat = uCovInvTL, dvec = rep( 0, nCoefTL ),
  Amat = t( restrTL$RMat ), bvec = - restrTL$RMat \%*\% uCoefTL + restrTL$rVec )
rCoefTL <- minDistTL$solution + uCoefTL
}

\keyword{models}
