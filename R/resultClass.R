###################################################################################
##' Constructor of Output class
##'
##' This class contains a result of a run.
##'
##' \describe{
##'   \item{proportion}{a vector of size K (number of clusters), containing the estimated mixture proportions.}
##'   \item{pi}{a matrix of size K*p (p= number of dimensions), containing the estimated probabilities of good paired comparison (dispersion parameters).}
##'   \item{mu}{a matrix with K rows containing the estimates modal rankings.}
##'   \item{ll}{the estimated log-likelihood.}
##'   \item{bic}{the estimated BIC criterion.}
##'   \item{icl}{the estimated ICL criterion.}
##'   \item{tik}{a matrix of size n*K (n= sample size) containing the posterior probability for the ith observation to belong to the kth cluster.}
##'   \item{partition}{a vector of size n containing the estimated partition of the observations.}
##'   \item{entropy}{a matrix of size n*2 containing the entropy of each individual and its estimated cluster.}
##'   \item{probability}{a matrix of size n*2 containing the probability of each individual and its estimated cluster.}
##'   \item{convergence}{if FALSE, no convergence has been achieved, no result is available.}
##'   \item{partial}{if FALSE, there is no partial rank in the data.}
##'   \item{partialRank}{a matrix containing the full rankings, observed or estimated (for partial rankings).}
##'   \item{distanceProp}{a list containing the distances between final estimation of the proportions and current estimation of the proportions at each iteration of the algorithm.}
##'   \item{distancePi}{a list containing the distances between final estimation of dispersion parameters and current estimation of the dispersion parameters at each iteration of the algorithm.}
##'   \item{distanceMu}{a list containing the distances between the final modal rankings estimation and current estimation of the modal rankings at each iteration of the algorithm.}
##'   \item{distanceZ}{a vector containing the Rand index between final partition and the current partition at each iteration of the algorithm.}
##'   \item{distancePartialRank}{a list containing the Kendall distances between final estimation of the partial rankings and current estimation of the partial rankings at each iteration of the algorithm.}
##'   \item{proportionInitial}{a vector containing the initialization of the proportions in the algorithm.}
##'   \item{piInitial}{a matrix containing the initialization of the probabilities of good paired comparison in the algorithm.}
##'   \item{muInitial}{a matrix containing the initialization of modal rankings in the algorithm.}
##'   \item{partialRankInitial}{a matrix containing the initialization of the partial rankings in the algorithm.}
##' }
##'
##'
##' @name Output-class
##' @rdname Output-class
## @exportClass Output
##'
setClass(
  Class="Output",
  representation=representation(
    proportion="numeric",
    pi="matrix",
    mu="matrix",
    ll="numeric",
    bic="numeric",
    icl="numeric",
    tik="matrix",
    partition="numeric",
	entropy="matrix",
	probability="matrix",
	convergence="logical",
	partial="logical",
	partialRank="matrix",
	distanceProp="list",
	distancePi="list",
	distanceMu="list",
	distanceZ="numeric",
	distancePartialRank="list",
	proportionInitial="numeric",
    piInitial="matrix",
    muInitial="matrix",
	partialRankInitial="matrix"
    ),
  prototype=prototype(
    proportion=numeric(0),
    pi=matrix(nrow=0,ncol=0),
    mu=matrix(nrow=0,ncol=0),
    ll=numeric(0),
    bic=numeric(0),
    icl=numeric(0),
    tik=matrix(nrow=0,ncol=0),
    partition=numeric(0),
	entropy=matrix(nrow=0,ncol=0),
	probability=matrix(nrow=0,ncol=0),
	convergence=logical(0),
	partial=logical(0),
	partialRank=matrix(nrow=0,ncol=0),
	distanceProp=list(),
	distancePi=list(),
	distanceMu=list(),
	distanceZ=numeric(0),
	distancePartialRank=list(),
	proportionInitial=numeric(0),
    piInitial=matrix(nrow=0,ncol=0),
    muInitial=matrix(nrow=0,ncol=0),
	partialRankInitial=matrix(nrow=0,ncol=0)
  )
)



###################################################################################
##' Constructor of Rankclust class
##'
##' This class contains results of rankclust function.
##'
##' \describe{
##'   \item{K}{vector containing the number of clusters.}
##'   \item{data}{data used in algorithm.}
##'   \item{criterion}{criterion defined in rankclust function to select the best result.}
##'   \item{convergence}{if 0, no convergence, no result available in results.}
##'   \item{results}{a list of the same length than K, containing Output objects.}
##' }
##'
##'
##' @name Rankclust-class
##' @rdname Rankclust-class
## @exportClass Rankclust
##'
setClass(
  Class="Rankclust",
  representation=representation(
    K="numeric",
    results="list",
	data="matrix",
    criterion="character",
	convergence="logical"
  ),
	prototype=prototype(
    results=list(),
	data=matrix(ncol=0,nrow=0),
    K=numeric(0),
    criterion="bic",
	convergence=logical(0)
  )

)


#' Getter method for rankclust output
#' 
#' This is overloading of square braces to extract values of various 
#' slots of the output from the function \code{\link{rankclust}}.
#' 
#' @param x object from which to extract element(s) or in which to replace element(s).
#' @param i the number of cluster of the element we want to extract.
#' 
#' @name [
#' 
#' @rdname getter-methods
#' @aliases [,Rankclust-method
#' 
setMethod(
  f="[",
  signature="Rankclust",
  definition=function(x,i,j,drop){
	if(x@convergence)
	{
		if(is.numeric(i))
		{
			if(i %in% x@K)
			{
				return(x@results[[which(x@K==i)]])       
			}
			else
				stop("Invalid number of cluster.") 
		}
		else
		{
			if(i=="bic")
			{
				bic=rep(NA,length(x@K))
				for(iter in 1:length(x@K))
				{
					if(x@results[[iter]]@convergence)
						bic[iter]=x@results[[iter]]@bic
				}
				return(bic)
			}
			else
			{
				if(i=="icl")
				{
					icl=rep(NA,length(x@K))
					for(iter in 1:length(x@K))
					{
						if(x@results[[iter]]@convergence)
							icl[iter]=x@results[[iter]]@icl
					}
					return(icl)
				}
				else
				{
					if(i=="ll")
					{
						ll=rep(NA,length(x@K))
						for(iter in 1:length(x@K))
						{
							if(x@results[[iter]]@convergence)
								ll[iter]=x@results[[iter]]@ll
						}
						return(ll)
					}
					else
					{
						stop("Invalid Name.")
					}
		
				}

			}
		}
    }
  }
)

#'
#' summary function.
#' 
#' This function This function gives the summary of output from \code{rankclust}.
#' 
#' @param object output object from \code{\link{rankclust}}.
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' @aliases summary summary,Rankclust-method
setMethod(
  f="summary",
  signature = "Rankclust",
  definition = function(object,...) {
	if(object@convergence)
	{
		if (object@criterion=="bic") 
		{
			BIC=c()
			for(i in object@K)
			{
				BIC=c(BIC,object@results[[which(object@K==i)]]@bic)
			}
			index=which(BIC==min(BIC))
		}
		else
		{
			ICL=c()
			for(i in object@K)
			{
				ICL=c(ICL,object@results[[which(object@K==i)]]@icl)
			}
			index=which(ICL==min(ICL))

		}

		cat("******************************************************************\n")
		cat("NUMBER OF CLUSTERS: ",object@K[index],"\n")
		if(object@criterion=="bic")  
			cat(object@criterion,"=",object[object@K[index]]@bic)
		else
			cat(object@criterion,"=",object[object@K[index]]@icl)
		cat("\nLoglikelihood =",object[object@K[index]]@ll)
		cat("\n\n*************************PARAMETERS*******************************\n")
		cat("Proportion:",object[object@K[index]]@proportion)
		cat("\nReference ranks mu:\n")
		print(object[object@K[index]]@mu)
		cat("\nProbabilities pi:\n")
		print(object[object@K[index]]@pi)
		cat("\n*************************CLUSTERING*******************************\n")
		cat("Ranks with the highest entropy for each cluster:\n")
		for(i in 1:object@K[index])
		{
			#classe=object[object@K[index]]@entropy[object[object@K[index]]@entropy[,2]==i,]
			classe=which(object[object@K[index]]@entropy[,2]==i)
			if(length(classe)!=0)
			{
				classe=classe[order(object[object@K[index]]@entropy[classe,1],decreasing=TRUE)][1:min(5,length(classe))]
				#if(object@algorithm=="SEM")
				print(cbind(object@data[classe,],object[object@K[index]]@entropy[classe,]))
				#else
				#	print(cbind(object@data[classe,-ncol(object@data)],object[object@K[index]]@entropy[classe,]))
				#if(length(classe)==2)
				#{
				#	best5=classe
				#	print(cbind(object@data[classe,-ncol(object@data)],best5[2:3]))
				#}
				#else
				#{
				#	best5=classe[order(classe[,1],decreasing=TRUE),][1:min(5,nrow(classe)),]
				#	print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
				#}
			}
	
		}
		#rm(best5)	
		cat("Ranks with the highest probability for each cluster:\n")
		for(i in 1:object@K[index])
		{
			#classe=object[object@K[index]]@probability[object[object@K[index]]@probability[,2]==i,]
			classe=which(object[object@K[index]]@probability[,2]==i)
			if(length(classe)!=0)
			{
				classe=classe[order(object[object@K[index]]@probability[classe,1],decreasing=TRUE)][1:min(5,length(classe))]
				#if(object@algorithm=="SEM")
				print(cbind(object@data[classe,],object[object@K[index]]@probability[classe,]))
				#else
				#	print(cbind(object@data[classe,-ncol(object@data)],object[object@K[index]]@probability[classe,]))
				#if(length(classe)==2)
				#{
				#	best5=classe
				#	print(cbind(object@data[best5[1],-ncol(object@data)],best5[2:3]))
				#}
				#else
				#{
				#	best5=classe[order(classe[,2],decreasing=TRUE),][1:min(5,nrow(classe)),]
				#	print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
				#}
				
			}	
		}   
	  	rm(classe)
		#rm(best5)
		if(object[object@K[index]]@partial)
		{
			cat("\n*************************PARTIAL RANK*****************************\n")
    		if(min(50,nrow(object[object@K[index]]@partialRank))==50)
      			cat("\nOnly the first 50 are printed, total length:",nrow(object[object@K[index]]@partialRank),"\n")
			print(object[object@K[index]]@partialRank[1:min(50,nrow(object[object@K[index]]@partialRank)),])
		}
		
		cat("\n******************************************************************\n")
	}
	else
		cat("\nNo convergence. Please retry\n")

  }  
)

#'
#' show function.
#' 
#' This function shows the elements of an object of class Output.
#' 
#' @param object of class Output.
#' 
#' @name show
#' @rdname show-methods
#' @docType methods
#' @exportMethod show
#' @aliases show show,Output-method
setMethod(
  f="show",
  signature = "Output",
  definition = function(object) {
    cat("ll=",object@ll)
    cat("\nbic =",object@bic)
    cat("\nicl =",object@icl)
    cat("\nproportion:",object@proportion)
    cat("\nmu:\n")
    print(object@mu)
    cat("\npi:\n")
    print(object@pi)
    cat("\npartition:\n")
    print(object@partition[1:min(50,length(object@partition))])
    if(min(50,length(object@partition))==50)
      cat("\nOnly the first 50 are printed, total length:",length(object@partition))
    cat("\ntik:\n")
    print(object@tik[1:min(50,nrow(object@tik)),])
    if(min(50,nrow(object@tik))==50)
      cat("\nOnly the first 50 rows are printed, total rows:",nrow(object@tik))
    }
)

#'
#' show function.
#' 
#' This function shows the output from \code{rankclust}.
#' 
#' @param object output object from \code{\link{rankclust}}.
#' 
#' @name show
#' @rdname show-methods
#' @docType methods
#' @exportMethod show
#' @aliases show show,Rankclust-method
#' 
setMethod(
  f="show",
  signature = "Rankclust",
  definition = function(object) {
    for(i in object@K)
	{
		cat("\n******************************************************************\n")
		cat("Number of clusters:",i)
		cat("\n******************************************************************\n")
		show(object@results[[which(object@K==i)]])
		cat("\n******************************************************************\n")
		
	}
    }
)
