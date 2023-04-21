\name{distRayEla}
\alias{distRayDeriv}
\alias{distRayEla}
\title{Elasticities and Derivatives of Ray-Based Input Distance Functions}

\description{
   Calculate distance elasticities and derivatives 
   based on ray-based input distance functions.
}

\usage{
distRayEla( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE, \dots )

distRayDeriv( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE, 
  numDeriv = FALSE, eps = 1e-6 ) 
}

\arguments{
  \item{xNames}{a vector of character strings 
    containing the names of the variables
    that indicate the input quantities.}
  \item{yNames}{a vector of two or more character strings 
    containing the names of the variable 
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
  \item{coef}{numeric vector containing the coefficients.}
  \item{data}{data frame containing the data.}
  \item{form}{a character string that indicates the functional form;
    currently, \code{"cd"} for the Cobb-Douglas functional form
    and \code{"tl"} for the Translog functional form are available.}
  \item{conDummy}{an optional numeric vector indicating the positions
    in argument \code{zNames} that are dummy variables
    so that quadratic terms of these variables are omitted.}
  \item{fixThetas}{logical value that indicates 
    whether undefined angles of the output  should be \sQuote{fixed}
    if the last two or more output quantities are zero 
    for some of the observations.}
  \item{numDeriv}{logical value that indicates 
    whether derivatives (and elasticities) should be calculated
    by numerical finite-difference differentiation.}
  \item{eps}{small positive value that is used as change in the variables
    when calculating derivatives (and elasticities)
    by numerical finite-difference differentiation.}
  \item{\dots}{further arguments of \code{distRayEla}
    are passed to \code{distRayDeriv}.}
}

\value{
   A list that will be described here later.
}

\author{Arne Henningsen and Juan José Price}

\examples{
# load and prepare data set 
data( appleProdFr86, package = "micEcon" )
appleProdFr86$qCap <- appleProdFr86$vCap / appleProdFr86$pCap
appleProdFr86$qLab <- appleProdFr86$vLab / appleProdFr86$pLab
appleProdFr86$qMat <- appleProdFr86$vMat / appleProdFr86$pMat

# estimate Translog ray-based input distance function
estTL <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86 )
summary( estTL )

# calculate elasticities
ela <- distRayEla( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estTL ), data = appleProdFr86 )
summary( ela )

# calculate derivatives
deriv <- distRayDeriv( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  coef = coef( estTL ), data = appleProdFr86 )
summary( deriv )
}

\keyword{models}
