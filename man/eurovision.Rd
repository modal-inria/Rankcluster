\docType{data}
\name{eurovision}
\alias{eurovision}
\title{Multidimensionnal partial rank data : eurovision}
\format{A list containing:
\describe{
\item{data}{ A matrix of size 34*48. Each row corresponds to the ordering representation of a multidimensionnal ranking.
Columns 1 to 8 correspond to the 2007 contest, columns 9 to 18 to the 2008 contest, etc...}

\item{frequency}{A matrix of size 34*49 containing the differents multidimensionnal rankings. The 48 first colmuns are the same as in data, and the last column contains the frequency (1 for all ranks).}

\item{m}{ a vector with the sizes of ranks for each dimension.}
}

}
\source{
  http://www.eurovision.tv
}
\description{
This dataset contains the ranking of the 8 common finalists of the Eurovision song contest from 2007 to 2012:

1: France, 2:Germany, 3:Greece, 4:Romania, 5:Russia, 6:Spain, 7:Ukrain, 8:United Kingdom.

The number of rankings is 33, corresponding to the 33 European countries having participated to this six editions of the contest.

All the rankings are partial since none country has ranked this 8 countries in its 10 preferences. Missing ranking elements are zeros.
	
}
\examples{
data(eurovision)
}
\keyword{datasets}

