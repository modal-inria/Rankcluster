\docType{data}
\name{words}
\alias{words}
\title{rank data : words}
\format{A list containing :
\describe{
\item{data}{A matrix of size 98*5 containing the 98 answers. Each row corresponding to one ranking (ordering representation).}
\item{frequency}{matrix of size 15*6. Each row corresponds to one of the 15 different observed rankings, and the last column contains the observation frequency.}
\item{m}{the size of the rankings (m=5).}
}
}
\source{
M.A. Fligner and J.S. Verducci. "Distance based ranking models". J. Roy. Statist. Soc. Ser. B, 48(3):359-369, 1986.
}
\description{
The data was collected under the auspices of the Graduate Record
Examination Board. A sample of 98 college students were asked to rank five words according to strength of association (least to most associated) with the target word "Idea":
1 = Thought, 2 = Play, 3 = Theory, 4 = Dream and 5 = Attention.
	
}
\examples{
data(words)
}
\keyword{datasets}

