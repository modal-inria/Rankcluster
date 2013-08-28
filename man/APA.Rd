\docType{data}
\name{APA}
\alias{APA}
\title{rank data : APA}
\format{A list containing :
\describe{
\item{data}{ matrix of size 5738x5 containing the 5738 observed full ranks in ordering representation.
The ordering representation x = (x_1,...,x_m) means that the candidate x_i is ranked in the ith position (i=1,...,m).
For example, if the ordering rank is (4,3,1,2,5), it means that judge ranks first the candidate 4, second the candidate 3, ...}

\item{frequency}{matrix of size 120x6. Each row corresponds to one of the different observed rank.
The first fifth columns contains the observed ranks (ordering representation) and the sixth column
contains the frequency of observation.}

\item{m}{ vector with the size of the ranks (5 here).}
}

}
\source{
  "Group representations in probability and statistics", P. Diaconis, 1988.
}
\description{
This dataset contains the 5738 full rankings resulting from the American Psychological Association (APA) presidential election of 1980.
For this election, members of APA had to rank five candidates in order of preference. 

For information, a total of 15449 votes have been registred for this election, but only the 5738 full rankings are reported in the APA dataset. Candidates 1 and 3 were research psychologists, candidates 4 and 5 were clinical psychologists and candidate 2 was a community psychologist.
	
}
\examples{
data(APA)
}
\keyword{datasets}

