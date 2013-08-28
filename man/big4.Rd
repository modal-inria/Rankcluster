\docType{data}
\name{big4}
\alias{big4}
\title{rank data : big4}
\format{A list containing :
\describe{
\item{data}{A matrix of size 21*8 containing the 21 Premier League seasons. Each row corresponding to one ranking (ordering representation).}
\item{frequency}{matrix of size 21*9. Each row corresponds to one of the 21 different observed rankings, and the last column contains the observation frequency.}
\item{m}{the size of the rankings (m=c(4,4) ).}
}
}
\source{

http://en.wikipedia.org/wiki/Big_Four_%28English_football%29#.22Big_Four.22_dominance_.282003.E2.80.9304_to_2008.E2.80.9309.29

http://www.uefa.com/memberassociations/uefarankings/club/index.html

}
\description{
This dataset is composed of the rankings (in ordering notation) of the "Big Four" English football teams (1: Manchester, 2: Liverpool, 3: Arsenal, 4: Chelsea) to the English Championship (Premier League) and according to the UEFA coefficients (statistics used in Europe for ranking and seeding teams in international competitions), from 1993 to 2003.
	
}
\examples{
data(big4)
}
\keyword{datasets}

