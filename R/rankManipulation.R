# convert a single rank
convertSingleRank <- function(x) {
  xb <- rep(0, length(x))
  ind <- x[x > 0]
  for (i in ind) {
    xb[i] <- which(x == i)
  }

  return(xb)
}

#' @title Change the representation of a rank
#'
#' @description convertRank converts a rank from its ranking representation to its ordering representation, and vice-versa.
#' The function does not work with partial ranking. The transformation to convert a rank from ordering to ranking
#' representation is the same that from ranking to ordering representation, there is no need to precise the representation
#' of rank x.
#'
#' @param x a rank (vector) datum either in its ranking or ordering representation.
#'
#' @return a rank (vector) in its ordering representation if its ranking representation has been given in input of
#' convertRank, and vice-versa.
#'
#'
#' @details The ranking representation r=(r_1,...,r_m) contains the ranks assigned to the objects,
#' and means that the ith object is in r_ith position.
#'
#' The ordering representation o=(o_1,...,o_m) means that object o_i is in the ith position.
#'
#' Let us consider the following example to illustrate both notations: a judge, which has to rank three holidays destinations
#' according to its preferences, O1 = Countryside, O2 =Mountain and O3 = Sea, ranks first Sea, second Countryside,
#' and last Mountain. The ordering result of the judge is o = (3, 1, 2) whereas the ranking result is r = (2, 3, 1).
#'
#'
#' @examples
#' x <- c(2, 3, 1, 4, 5)
#' convertRank(x)
#'
#' @author Julien Jacques
#'
#' @export
convertRank <- function(x) {
  if (is.matrix(x)) {
    return(t(apply(x, 1, convertSingleRank)))
  } else {
    return(convertSingleRank(x))
  }
}

# '
# ' This function checks if data are correct.
# '
# '
# ' @title Check the data
# ' @param X a matrix containing ranks
# ' @param m a vector composed of the sizes of the rankings of each dimension
# ' (default value is the number of column of the matrix data).
# '
# ' @return a list containing for each dimension, numbers of rows with problem.
# '
# ' @examples
# ' data(big4)
# ' #add a bad rank
# ' big4$data[1,1:4] = c(1,5,2,3)
# '
# ' res=checkData(big4$data,big4$m)
# ' print(res)
# '
# ' @export
# checkData = function(X,m=length(X))
# {
#   if(!is.matrix(X))
#   {
#     X=as.matrix(X)
#   }
#   if(sum(m)!=ncol(X))
#     stop("the number of columns of X does not match with m.")
#
#   d=length(m)
#   n=nrow(X)
#   cm=cumsum(c(0,m))
#   pb=list()
#   for(i in 1:d)
#   {
#     check=apply(X[,(1+cm[i]):cm[i+1],drop=FALSE],1,checkTiePartialRank,m[i])
#     if(sum(check)!=n)
#     {
#       indfalse=which(check==0)
#       pb[[i]]=indfalse
#     }
#   }
#
#   return(pb)
# }

# checkRank  check if a vector is a rank
checkRank <- function(x, m = length(x)) {
  return(sum(sort(x) == (1:m)) == m)
}

# checkPartialRank check if a vector is a partial rank
# missing element : 0
checkPartialRank <- function(x, m = length(x)) {
  return((length(x[x <= m]) == m) && (length(x[x >= 0]) == m) && (length(unique(x[x != 0])) == length(x[x != 0])))
}


# checkPartialRank check if a vector is a partial rank
checkTiePartialRank <- function(x, m = length(x)) {
  return((length(x[x <= m]) == m) && (length(x[x >= 0]) == m))
}

# completeRank complete partial that have only one missing element
completeRank <- function(x) {
  if (length(x[x == 0]) == 1) {
    m <- length(x)
    a <- 1:m
    a[x[x != 0]] <- 0
    x[x == 0] <- a[a != 0]
  }
  return(x)
}


#' @title Convert data storage
#'
#' @description This function takes in input a matrix containing all the observed ranks (a rank can be repeated)
#' and returns a matrix containing all the different observed ranks with their observation frequencies (in the last column).
#'
#' @param X a matrix containing ranks.
#' @param m a vector with the size of ranks of each dimension.
#'
#' @return A matrix containing each different observed ranks with its observation frequencies in the last column.
#'
#' @examples
#' X <- matrix(1:4, ncol = 4, nrow = 5, byrow = TRUE)
#' Y <- frequence(X)
#' Y
#'
#' @author Quentin Grimonprez
#'
#' @seealso \link{unfrequence}
#'
#' @export
frequence <- function(X, m = ncol(X)) {
  if (missing(X)) {
    stop("X is missing")
  }
  if (!is.numeric(X) || !is.matrix(X)) {
    stop("X must be a matrix of positive integer")
  }
  if (length(X[X >= 0]) != length(X)) {
    stop("X must be a matrix of positive integer")
  }
  if (!is.vector(m, mode = "numeric")) {
    stop("m must be a (vector of) integer strictly greater than 1")
  }
  if (length(m) != length(m[m > 1])) {
    stop("m must be a (vector of) integer strictly greater than 1")
  }

  if (length(m) == 1) {
    if (m != ncol(X)) {
      print(paste0(
        "You put m=", m, ", but X has ", ncol(X), " columns(rank of size ", ncol(X) - 1, " and 1 for the frequence)."
      ))
      print(paste0("The algorithm will continue with m=", ncol(X) - 1))
    }
  }

  res <- .Call("freqMultiR", X, m, PACKAGE = "Rankcluster")

  data <- matrix(0, ncol = length(res$data[[1]]) + 1, nrow = length(res$data))
  for (i in seq_len(nrow(data))) {
    data[i, ] <- c(res$data[[i]], res$freq[[i]])
  }


  return(data)
}


#' @title Convert data
#'
#' @description This function takes in input a matrix in which the m first columns are the different observed ranks and
#' the last column contains the observation frequency, and returns a matrix containing all the ranks (ranks with frequency>1
#' are repeated).
#'
#' @param data a matrix containing rankings and observation frequency.
#'
#' @return a matrix containing all the rankings.
#'
#' @examples
#' data(quiz)
#' Y <- unfrequence(quiz$frequency)
#' Y
#'
#' @seealso \link{frequence}
#'
#' @export
unfrequence <- function(data) {
  X <- matrix(ncol = ncol(data) - 1, nrow = sum(data[, ncol(data)]))
  colnames(X) <- colnames(data)[-ncol(data)]
  compteur <- 1
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(data[i, ncol(data)])) {
      X[compteur, ] <- data[i, -ncol(data)]
      compteur <- compteur + 1
    }
  }

  return(X)
}
