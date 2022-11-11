#' @title Simulate a sample of ISR(pi,mu)
#'
#' @description This function simulates univariate rankings data (ordering representation) according to the ISR(pi,mu).
#'
#' @param n size of the sample.
#' @param pi dispersion parameter: probability of correct paired comparison according to mu.
#' @param mu position parameter: modal ranking in ordering representation.
#' @return a matrix with simulated ranks.
#'
#' @details
#' The ranking representation r=(r_1,...,r_m) contains the
#' ranks assigned to the objects, and means that the ith
#' object is in r_ith position.
#'
#' The ordering representation o=(o_1,...,o_m) means that
#' object o_i is in the ith position.
#'
#' Let us consider the following example to illustrate both
#' notations: a judge, which has to rank three holidays
#' destinations according to its preferences, O1 =
#'   Countryside, O2 =Mountain and O3 = Sea, ranks first Sea,
#' second Countryside, and last Mountain. The ordering
#' result of the judge is o = (3, 1, 2) whereas the ranking
#' result is r = (2, 3, 1).
#'
#' You can see the \link{convertRank} function to convert the simulated ranking from ordering to ranking representation.
#'
#' @references
#' [1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm,
#' Computational Statistics and Data Analysis, 58, 162-176.
#'
#' @examples
#' x <- simulISR(30, 0.8, 1:4)
#'
#' @author Julien Jacques
#'
#' @export
simulISR <- function(n, pi, mu) {
  if (missing(n)) {
    stop("n is missing")
  }
  if (missing(mu)) {
    stop("mu is missing")
  }
  if (missing(pi)) {
    stop("pi is missing")
  }

  if (!is.numeric(n) || (length(n) > 1)) {
    stop("n must be a strictly positive integer")
  }
  if ((n != round(n)) || (n <= 0)) {
    stop("n must be a strictly positive integer")
  }

  if (!is.numeric(pi) || (length(pi) > 1)) {
    stop("pi must be a real between 0 and 1")
  }
  if ((pi > 1) || (pi < 0)) {
    stop("pi must be a real between 0 and 1")
  }

  if (!is.vector(mu, mode = "numeric")) {
    stop("mu must be a complete rank")
  }
  if (!checkRank(mu)) {
    stop("mu must be a complete rank")
  }



  res <- .Call("simulISRR", n, length(mu), mu, pi, PACKAGE = "Rankcluster")

  return(res)
}


#' @title Probability computation
#'
#' @description It computes the probability of a (multivariate) rank x according to a ISR(mu, pi).
#'
#' @param x a vector or a matrix containing the rankings in ranking notation (see Details or \link{convertRank} function).
#' The rankings of each dimension are placed end to end. \code{x} must contain only full ranking (no partial or tied).
#' @param pi a vector of size \code{p=length(m)}, where \code{p} is the number of dimension, containing the probabilities of
#' a good comparison of the model (dispersion parameters).
#' @param mu a vector of length \code{sum(m)} containing the modal ranks in ranking notation (see Details or
#' \link{convertRank} function).
#' The rankings of each dimension are placed end to end. \code{mu} must contain only full ranking (no partial or tied).
#' @param m a vector containing the size of ranks for each dimension.
#'
#' @return the probability of \code{x} according to a ISR(mu, pi).
#'
#' @details
#' The ranks have to be given to the package in the ranking notation (see \link{convertRank} function),
#' with the following convention:
#'
#' - missing positions are replaced by 0
#'
#' - tied are replaced by the lowest position they share"
#'
#'
#'   The ranking representation r=(r_1,...,r_m) contains the
#' ranks assigned to the objects, and means that the ith
#' object is in r_ith position.
#'
#' The ordering representation o=(o_1,...,o_m) means that object
#' o_i is in the ith position.
#'
#' Let us consider the following example to illustrate both
#' notations: a judge, which has to rank three holidays
#' destinations according to its preferences, O1 =
#'   Countryside, O2 =Mountain and O3 = Sea, ranks first Sea,
#' second Countryside, and last Mountain. The ordering
#' result of the judge is o = (3, 1, 2) whereas the ranking
#' result is r = (2, 3, 1).
#'
#'
#' @examples
#' m <- c(4, 5)
#' x <- mu <- matrix(nrow = 1, ncol = 9)
#' x[1:4] <- c(1, 4, 2, 3)
#' x[5:9] <- c(3, 5, 2, 4, 1)
#' mu[1:4] <- 1:4
#' mu[5:9] <- c(3, 5, 4, 2, 1)
#' pi <- c(0.75, 0.82)
#'
#' prob <- probability(x, mu, pi, m)
#' prob
#'
#' @author Quentin Grimonprez
#'
#' @export
#'
probability <- function(x, mu, pi, m = length(mu)) {
  x <- matrix(x, ncol = sum(m))
  mu <- matrix(mu, ncol = sum(m))

  checkProbability(x, mu, pi, m)

  # convert to ordering
  for (i in seq_along(m)) {
    mu[1, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])] <- convertRank(
      mu[1, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])]
    )
    for (j in seq_len(nrow(x))) {
      x[j, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])] <- convertRank(
        x[j, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])]
      )
    }
  }

  proba <- rep(NA, nrow(x))
  for (j in seq_along(proba)) {
    proba[j] <- .Call("computeProba", x[j, , drop = FALSE], mu, pi, m, PACKAGE = "Rankcluster")
  }

  return(proba)
}


checkProbability <- function(x, mu, pi, m) {
  ### check parameters
  if (missing(x)) {
    stop("x is missing.")
  }
  if (missing(mu)) {
    stop("mu is missing.")
  }
  if (missing(pi)) {
    stop("pi is missing.")
  }

  # x
  if (!(is.vector(x) || is.matrix(x))) {
    stop("x must be either a matrix or a vector.")
  }
  if (is.vector(x)) {
    x <- t(as.matrix(x))
  }
  if (!is.numeric(x)) {
    stop("x must be either a matrix or a vector of integer.")
  }

  # mu
  if (!(is.vector(mu) || is.matrix(mu))) {
    stop("mu must be either a matrix or a vector.")
  }
  if (is.vector(mu)) {
    mu <- t(as.matrix(mu))
  }
  if (!is.numeric(mu)) {
    stop("mu must be either a matrix or a vector of integer.")
  }

  # pi
  if (!is.numeric(pi)) {
    stop("pi must be a vector of probabilities.")
  }
  if (!is.vector(pi)) {
    stop("pi must be a vector of probabilities.")
  }
  if ((min(pi) < 0) || max(pi) > 1) {
    stop("pi must be a vector of probabilities.")
  }

  # m
  if (!is.numeric(m)) {
    stop("m must be a vector of integer.")
  }
  if (!is.vector(m)) {
    stop("m must be a vector of integer.")
  }
  if (sum(unlist(lapply(m, is.wholenumber))) != length(m)) {
    stop("m contains non integer.")
  }
  if (sum(m) != ncol(x)) {
    stop("sum(m) and the length of x do not match.")
  }
  if (sum(m) != length(mu)) {
    stop("sum(m) and the length of mu do not match.")
  }
  if (length(m) != length(pi)) {
    stop("the length of pi and m do not match.")
  }

  # check if mu contains ranks
  for (i in seq_along(m)) {
    if (!checkRank(mu[, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])], m[i])) {
      stop("mu is not correct.")
    }
    for (j in seq_len(nrow(x))) {
      if (!checkRank(x[j, (1 + cumsum(c(0, m))[i]):(cumsum(c(0, m))[i + 1])], m[i])) {
        stop("x is not correct.")
      }
    }
  }
}