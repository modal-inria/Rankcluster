\name{simulISR}
\alias{simulISR}
\title{simulate a sample of ISR(pi,mu)}
\usage{
  simulISR(n, pi, mu)
}
\arguments{
  \item{n}{size of the sample.}

  \item{pi}{dispersion parameter: probability of correct
  paired comparaison according to mu.}

  \item{mu}{position parameter: modal ranking in ordering
  representation.}
}
\value{
  a matrix with simulated ranks.
}
\description{
  This function simulates univariate rankings data
  according to the ISR(pi,mu).
}
\examples{
x=simulISR(30,0.8,1:4)
}
\author{
  Julien Jacques
}
\references{
  [1] C.Biernacki and J.Jacques (2013), A generative model
  for rank data based on sorting algorithm, Computational
  Statistics and Data Analysis, 58, 162-176.
}

