\name{distRayEst}
\alias{distRayEst}
\title{Estimate a Ray-Based Input Distance Function}

\description{
   Empirically analyse multiple-input-multiple-output production technologies 
   by estimating a ray-based input distance function.
}

\usage{
distRayEst( xNames, yNames, zNames = NULL, sNames = NULL, 
  eNames = NULL, data, form = "tl", method = "sfa", fixThetas = FALSE, \dots )
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
  \item{eNames}{an optional vector of character strings 
    containing ???.}
  \item{data}{data frame containing the data.}
  \item{form}{a character string that indicates the functional form;
    currently, \code{"cd"} for the Cobb-Douglas functional form
    and \code{"tl"} for the Translog functional form are available.}
  \item{method}{a character string that indicates the estimation method;
    currently, \code{"lm"} for the Ordinary Least-Squares method
    (using the \code{lm()} function of the \code{stats} package)
    and \code{"sfa"} for stochastic frontier estimation
    (using the \code{sfacross()} function of the \code{sfaR} package)
    are available.}
  \item{fixThetas}{logical value that indicates 
    whether undefined angles of the output  should be \sQuote{fixed}
    if the last two or more output quantities are zero 
    for the some of the observations.}
  \item{\dots}{further arguments of \code{distRayEst}
    are passed to \code{\link[stats]{lm}}
    or \code{\link[sfaR]{sfacross}}.}
}

\value{
   A list that will be described here later.
}

\author{Arne Henningsen and Juan José Price}

\examples{
   data( germanFarms, package = "micEcon" )
   
   # to be added later
}

\keyword{models}
