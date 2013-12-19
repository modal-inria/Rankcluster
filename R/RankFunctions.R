#'convertRank converts a rank from its ranking representation to its ordering representation, and vice-versa. The function does not work with partial ranking.
#'The transformation to convert a rank from ordering to ranking representation is the same that from ranking to ordering representation, there is no need to precise the representation of rank x.
#'
#'The ranking representation r=(r_1,...,r_m) contains the ranks assigned to the objects,
#'and means that the ith object is in r_ith position.
#'
#'The ordering representation o=(o_1,...,o_m) means that object o_i is in the ith position. 
#'
#'
#'Let us consider the following example to illustrate both notations: a judge, which has to rank three holidays destinations according to its preferences, O1 = Countryside, O2 =Mountain and O3 = Sea, ranks first Sea, second Countryside, and last Mountain. The ordering result of the judge is o = (3, 1, 2) whereas the ranking result is r = (2, 3, 1).
#' @useDynLib Rankcluster
#' @title change the representation of a rank
#' @author Julien Jacques
#' @param x a rank (vector) datum either in its ranking or ordering representation.
#' @return a rank (vector) in its ordering representation if its ranking representation has been given in input of convertRank, and vice-versa.
#' @examples
#' x=c(2,3,1,4,5)
#' convertRank(x)
#' @export

convertRank <- function(x){
	if(is.matrix(x))
	{
		return(t(apply(x,1,FUN=function(x) {sort.int(x,index.return=1)$ix} ) ) )
	}
	else
	{
		return(sort.int(x,index.return=1)$ix)
	}	
}

# checkRank  check if a vector is a rank


checkRank <- function(x,m=length(x))
{
	if(sum(sort(x)==(1:m))==m)
		return(TRUE)
	else
		return(FALSE)	
}

# checkPartialRank check if a vector is a rank


checkPartialRank <- function(x,m=length(x))
{
	if((length(x[x<=m])==m)&& (length(x[x>=0])==m) && (length(unique(x[x!=0]))==length(x[x!=0])))
		return(TRUE)
	else
		return(FALSE)
}

# completeRank complete partial that have only one missing element


completeRank <-function(x)
{
	if(length(x[x==0])==1)
	{	
		m=length(x)
		a=1:m
		a[x[x!=0]]=0
		x[x==0]=a[a!=0]
	}
	return(x)
}

#' This function takes in input a matrix containing all the observed ranks (a rank can be repeated) and returns a matrix containing all the different observed ranks with their observation frequencies (in the last column).
#' @title Convert data storage
#' @author Quentin Grimonprez
#' @param X a matrix containing ranks.
#' @param m a vector with the size of ranks of each dimension.
#' @return A matrix containing each different observed ranks with its observation frequencies in the last column.
#' @examples
#' X=matrix(1:4,ncol=4,nrow=5,byrow=TRUE)
#' Y=frequence(X)
#' Y
#' @export
frequence <-function(X,m=ncol(X))
{
	if(missing(X))
		stop("X is missing")
	if(!is.numeric(X) || !is.matrix(X))
		stop("X must be a matrix of positive integer")
	if(length(X[X>=0])!=length(X))
		stop("X must be a matrix of positive integer")
	if(!is.vector(m,mode="numeric"))
		stop("m must be a (vector of) integer strictly greater than 1")
	if(length(m)!=length(m[m>1]))
		stop("m must be a (vector of) integer strictly greater than 1")

	if(length(m)==1)
	{
		if(m!=ncol(X))
		{
			print(paste0("You put m=",m,", but X has ",ncol(X)," columns(rank of size ",ncol(X)-1," and 1 for the frequence)."))
 			print(paste0("The algorithm will continue with m=",ncol(X)-1))
		}
	}

	res=.Call("freqMultiR",X,m,PACKAGE="Rankcluster")
	
	data=matrix(0,ncol=length(res$data[[1]])+1,nrow=length(res$data))
	for(i in 1:nrow(data))
		data[i,]=c(res$data[[i]],res$freq[[i]])
	

	return(data)

}

#' This function simulates univariate rankings data according to the ISR(pi,mu).
#' @title simulate a sample of ISR(pi,mu)
#' @author Julien Jacques
#' @param n size of the sample.
#' @param pi dispersion parameter: probability of correct paired comparaison according to mu.
#' @param mu position parameter: modal ranking in ordering representation.
#' @return a matrix with simulated ranks.
#' @references 
#' [1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176.
#' @examples
#' x=simulISR(30,0.8,1:4)
#' @export
simulISR <-function(n,pi,mu)
{
	if(missing(n))
		stop("n is missing")
	if(missing(mu))
		stop("mu is missing")
	if(missing(pi))
		stop("pi is missing")

	if(!is.numeric(n) || (length(n)>1))
		stop("n must be a strictly positive integer")
	if( (n!=round(n)) || (n<=0))
		stop("n must be a strictly positive integer")

	if(!is.numeric(pi) || (length(pi)>1))
		stop("pi must be a real between 0 and 1")
	if( (pi>1) || (pi<0))
		stop("pi must be a real between 0 and 1")

	if(!is.vector(mu,mode="numeric"))
		stop("mu must be a complete rank")
	if(!checkRank(mu))
		stop("mu must be a complete rank")



	res=.Call("simulISRR",n,length(mu),mu,pi,PACKAGE="Rankcluster")

	return(res)
}

#' This function takes in input a matrix in which the m first columns are the different observed ranks and the last column contains the observation frequency, and returns a matrix containing all the ranks (ranks with frequency>1 are repeated).
#' @title Convert data
#' @param data a matrix containing rankings and observation frequency.
#' @return a matrix containing all the rankings.
#' @examples
#' data(quiz)
#' Y=unfrequence(quiz$frequency)
#' Y
#' @export
unfrequence=function(data)
{
  X=matrix(ncol=ncol(data)-1,nrow=sum(data[,ncol(data)]))
  colnames(X)=colnames(data)[-ncol(data)]
  compteur=1
  for(i in 1:nrow(data))
    for(j in 1:data[i,ncol(data)])
    {
      X[compteur,]=data[i,-ncol(data)]
      compteur=compteur+1
    }
  return(X)    
}
