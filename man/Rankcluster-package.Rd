\name{Rankcluster-package}
\alias{Rankcluster-package}
\alias{Rankcluster}
\docType{package}
\title{
Model-based clustering for multivariate partial ranking data
}
\description{

This package proposes a model-based clustering algorithm for ranking data. 
Multivariate rankings as well as partial rankings are taken into account.
This algorithm is based on an extension of the Insertion Sorting Rank (ISR) model for ranking data, which is a meaningful
and effective model parametrized by a position parameter (the modal ranking, quoted by mu) and a dispersion parameter (quoted by pi).
The heterogeneity of the rank population is modelled by a mixture of ISR, whereas conditional independence assumption is considered for multivariate rankings.

}


\details{
\tabular{ll}{
Package: \tab Rankcluster\cr
Type: \tab Package\cr
Version: \tab 0.91.5\cr
Date: \tab 2013-11-03\cr
License: \tab GPL (>=2) \cr
}
}

\author{
Maintainer: Quentin Grimonprez <quentin.grimonprez@inria.fr>
}

\references{
[1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176.


[2] J.Jacques and C.Biernacki (2012), Model-based clustering for multivariate partial ranking data, Inria Research Report n 8113.
}

\keyword{ package }

