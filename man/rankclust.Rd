\name{rankclust}
\alias{rankclust}
\title{model-based clustering for multivariate partial ranking}
\usage{
  rankclust(data, m = ncol(data), K = 1, criterion = "bic",
    Qsem = 100, Bsem = 20, RjSE = m * (m - 1)/2,
    RjM = m * (m - 1)/2, Ql = 500, Bl = 100, maxTry = 3,
    run = 1, detail = FALSE)
}
\arguments{
  \item{data}{a matrix in which each row is a ranking
  (partial or not; for partial ranking, missing elements
  must be 0). For multivariate rankings, the rankings of
  each dimension are placed end to end in each row.}

  \item{m}{a vector composed of the sizes of the rankings
  of each dimension (default value is the number of column
  of the matrix data).}

  \item{K}{an integer or a vector of integer with the
  number of clusters.}

  \item{criterion}{criterion "bic" or "icl", criterion to
  minimize for selecting the number of clusters.}

  \item{Qsem}{the total number of iterations for the SEM
  algorithm (defaut value=40).}

  \item{Bsem}{burn-in period for SEM algorithm (default
  value=10).}

  \item{RjSE}{a vector containing, for each dimension, the
  number of iterations of the Gibbs sampler used both in
  the SE step for partial rankings and for the presentation
  orders generation (default value=mj(mj-1)/2).}

  \item{RjM}{a vector containing, for each dimension, the
  number of iterations of the Gibbs sampler used in the M
  step (default value=mj(mj-1)/2)}

  \item{Ql}{number of iterations of the Gibbs sampler for
  estimation of log-likelihood (default value=100).}

  \item{Bl}{burn-in period for estimation of log-likelihood
  (default value=50).}

  \item{maxTry}{maximum number of restarts of the SEM-Gibbs
  algorithm in the case of non convergence (default
  value=3).}

  \item{run}{number of runs of the algorithm for each value
  of K.}

  \item{detail}{boolean, if TRUE, time and others
  informations will be print during the process (default
  value FALSE).}
}
\value{
  An object of class rankclust.

  For example : res=rankclust(data,K=1:2,m=m)

  You can access the result by res[number of
  groups]@slotName where slotName is an element of the
  class Output.
}
\description{
  This functions estimates a clustering of ranking data,
  potentially multivariate and partial, based on a mixture
  of multivariate ISR model [2]. By specifying only one
  cluster, the function performs a modelling of the ranking
  data using the multivariate ISR model. The estimation is
  performed thanks to a SEM-Gibbs algorithm in the general
  case.
}
\examples{
data(big4)
result=rankclust(big4$data,K=2,m=big4$m,Ql=200,Bl=100,maxTry=2)
}
\author{
  Quentin Grimonprez
}
\references{
  [1] C.Biernacki and J.Jacques (2013), A generative model
  for rank data based on sorting algorithm, Computational
  Statistics and Data Analysis, 58, 162-176.

  [2] J.Jacques and C.Biernacki (2012), Model-based
  clustering for multivariate partial ranking data, Inria
  Research Report n 8113.
}
\seealso{
  See \code{\link{Output-class}} and
  \code{\link{Rankclust-class}} for available output.
}

