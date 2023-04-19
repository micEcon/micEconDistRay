\name{MuseumsDk}
\docType{data}
\alias{MuseumsDk}

\title{Data on Museums in Denmark}

\description{
  The \code{MuseumsDk} data set is a balanced panel data set 
  of 93 state-recognized museums in Denmark over a six years
  (2012 and 2014-2018; 2013 is unavailable).
}

\usage{data( "MuseumsDk" )}

\format{
  This data frame contains the following columns:
    \describe{
      \item{museum}{Name of the museum.}
      \item{type}{Type of museum 
        (Kulturhistorisk museum = cultural history museum;
        Kunstmuseer = arts museum;
        Naturhistorisk museum = natural history museum;
        Blandet museum = mixed museum).}
      \item{munic}{Municipality, in which the museum is located.}
      \item{yr}{Year of the observation.}
      \item{units}{Number of visit sites.}
      \item{resp}{whether or not the museum has special responsibilities
        (0 = no special responsibiliries; 1 = at least one special responsiility).}
      \item{vis}{Number of (physical) visitors.}
      \item{aarc}{Number of articles published (archeology).}
      \item{ach}{Number of articles published (cultural history).}
      \item{aah}{Number of articles published (art history).}
      \item{anh}{Number of articles published (natural history).}
      \item{exh}{Number of temporary exhibitions.}
      \item{edu}{Number of primary school classes on educational visits to the museum	.}
      \item{ev}{Number of events other than exhibitions.}
      \item{ftesc}{Scientific labor (full-time equivalents).}
      \item{ftensc}{Non-scientific labor (full-time equivalents).}
      \item{expProperty}{Running and maintenance costs [1,000 DKK].}
      \item{expCons}{conversation expenditure [1,000 DKK].}
      \item{ipc}{Consumer price index in Denmark 
        (the value for year 2014 is set to 1).}
    }
}

\source{
  A subset of this data set is used for the empirical analysis in
  Price \& Henningsen (forthcoming).
  It has been obtained from Statistics Denmark.
}

\references{
  Price, J.J. \& Henningsen, A. (forthcoming):
  A Ray-Based Input Distance Function to Model Zero-Valued Output Quantities: 
  Derivation and an Empirical Application.
  Journal of Productivity Analysis.
}

\keyword{datasets}
