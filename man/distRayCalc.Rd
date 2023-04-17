\name{distRayCalc}
\alias{distRayCalc}
\title{Calculate the Dependent Variable of a Ray-Based Input Distance Function}

\description{
   Calculate the dependent variable (logarithmic distance) 
   of a ray-based input distance function.
}

\usage{
distRayCalc( xNames, yNames, zNames = NULL, sNames = NULL, 
  coef, data, form = "tl", conDummy = NULL, fixThetas = FALSE )
}

\arguments{
  \item{xNames}{a vector of character strings 
    containing the names of the variabels
    that indicate the input quantities.}
  \item{yNames}{a vector of two or more character strings 
    containing the names of the variable 
    that indicate the output quantities.}
  \item{zNames}{an optional vector of character strings 
    containing the names of \sQuote{environmental} variabes,
    i.e., variables that affect the production possibility set
    (i.e., the feasible combinations of input-output quantities)
    that---in case of a Translog functional form---should be interacted 
    with the input quantities and the angels of the output vector.}
  \item{sNames}{an optional vector of character strings 
    containing the names of \sQuote{environmental} variabes,
    i.e., variables that affect the production possibility set
    (i.e., the feasible combinations of input-output quantities)
    that---in case of a Translog functional form---should not be interacted 
    with the input quantities and the angels of the output vector.}
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
    for the some of the observations.}
}

\value{
   A vector.
}

\author{Arne Henningsen and Juan José Price}

\examples{
# load and prepare data set 
data( appleProdFr86, package = "micEcon" )
appleProdFr86$qCap <- appleProdFr86$vCap / appleProdFr86$pCap
appleProdFr86$qLab <- appleProdFr86$vLab / appleProdFr86$pLab
appleProdFr86$qMat <- appleProdFr86$vMat / appleProdFr86$pMat

# Cobb-Douglas ray-based input distance function (with manually set parameters)
appleProdFr86$logDistCD <- distRayCalc( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ), data = appleProdFr86,
  coef = c( "(Intercept)" = -11.116, alpha_1 = 0.082, alpha_2 = 0.615,
  beta_1 = -0.031, beta_2 = -0.388 ), form = "cd" )
summary( appleProdFr86$logDistCD )

# Translog ray-based input distance function (with estimated parameters)
estTL <- distRayEst( xNames = c( "qCap", "qLab", "qMat" ),
  yNames = c( "qApples", "qOtherOut" ),
  data = appleProdFr86 )
appleProdFr86$logDistTL <- distRayCalc( xNames = c( "qCap", "qLab", "qMat" ),
  yNames= c( "qApples", "qOtherOut" ), data = appleProdFr86,
  coef = coef( estTL ) )
summary( appleProdFr86$logDistTL )
}

\keyword{models}
