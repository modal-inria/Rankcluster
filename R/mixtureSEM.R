
# mixtureSEM
# @title algorithm SEM for 
# @author Grimonprez Quentin
# @param X matrix where each row is a rank and the last column contains the frequencies
# @param g number of groups
# @param m a vector with the size of ranks for each dimension
# @param maxIt the maximum number of iteration of the algorithm

# @param Qsem the total number of iterations for the SEM algorithm (defaut value=40)
# @param Bsem burn-in period for SEM algorithm (default value=10)
# @param RjSE a vector containing the number of iteration for each dimension of the Gibbs algorithm in the SE step for generate partial ranks and orders of presentation(only for SEM algorithm, default value=m(m-1)/2)
# @param RjM a vector containing the number of iterations for each dimension for the Gibbs Sampler in the M step(only for SEM algorithm, default value=m(m-1)/2)
# @param Ql number of iterations of the Gibbs sampler for estimation of log-likelihood (only for SEM algorithm, default value=100)
# @param Bl burn-in period for estimation of log-likelihood (only for SEM algorithm, default value=50)
# @param detail boolean, if TRUE, time and others informations will be print during the process (default value FALSE)

# @return an object containing the refererence rank mu, the probability pi of a correct comparaison, proportion, conditionnal probability of belonging to each cluster (tik), the loglikelihood, the partition, the BIC and the ICL
# @references "Model-based clustering for multivariate partial ranking data", J. Jacques, C. Biernacki
# @examples
# data(APA)#m=5
# mixtureSEM(APA,2)
# @export

mixtureSEM<-function(X,g,m,Qsem,Bsem,Ql,Bl,RjSE,RjM,maxTry,run,detail)
{
	
	n=nrow(X)
	d=length(m)
	if(ncol(X)!=sum(m))
		stop(paste0("the number of column of X (",ncol(X),") does not match to the sum of vector m (",sum(m),")."))


	#Verification des donnees
	for(i in 1:d)
	{
		if(sum(apply(X[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkPartialRank,m[i]))!=n)
			stop("Data are not correct")
	}

	res=.Call("semR",X,m,g,Qsem,Bsem,Ql,Bl,RjSE,RjM,maxTry,run,detail,PACKAGE="Rankcluster")
	#récupération des résultats
	if(res$stock[1]==1)#si convergence
	{

		res$referenceRank=tliste3d2mat(res$referenceRank)
		res$initMu=tliste3d2mat(res$initMu)
		res$p=liste2d2matgd(res$p)
		res$initPi=liste2d2matgd(res$initPi)
		res$cluster=res$cluster+1

		res$entropy=cbind(res$entropy,res$cluster)
		res$probability=cbind(res$probability,res$cluster)
		colnames(res$entropy)=c("entropy","cluster")		
		colnames(res$probability)=c("probability","cluster")

		res$distMu=liste3d2listematgd(res$distMu)
		res$distP=liste3d2listematgd(res$distP)


		if(res$stock[2]==1)#si il y a des données partielles
		{
			#print(res$partialRank)#liste[[dim]][[ind]] vector
			#print(res$initPartialRank)#liste[[dim]][[ind]] vector
			#print(res$distPartialRank[[1]])#liste 30 iter
			#print(res$indexPartialData)

			res$partialRank=tliste3d2mat(res$partialRank)##proba a rajouté 

			rownames(res$partialRank)=rep("",nrow(res$partialRank))#enlever les cl1...
			#colnames(res$partialRank)[1]="Index"
			#colnames(res$partialRank)[ncol(res$rangPartial)]="Probability"

			res$initPartialRank=tliste3d2mat(res$initPartialRank)
			#colnames(res$initPartialRank)[1]="Index"
			rownames(res$initPartialRank)=rep("",nrow(res$initPartialRank))

			res$distPartialRank=lapply(res$distPartialRank,FUN=function(x){listedistPartiel(x)})
			result=new(Class="Output",
				bic=res$stock[4],
				icl=res$stock[5],
				ll=res$stock[3],
				proportion=res$proportion,
				pi=res$p,
				mu=res$referenceRank,
				tik=res$tik,
				partition=res$cluster,
				entropy=res$entropy,
				probability=res$probability,
				convergence=TRUE,
				partial=TRUE,
				partialRank=res$partialRank,
				distanceZ=res$distZ,
				distanceMu=res$distMu,
				distanceProp=res$distProp,
				distancePi=res$distP,
				distancePartialRank=res$distPartialRank,
				piInitial=res$initPi,
   				muInitial=res$initMu,
				partialRankInitial=res$initPartialRank,
				proportionInitial=res$initProportion)
		}
		else
		{
			
			result=new(Class="Output",
				bic=res$stock[4],
				icl=res$stock[5],
				ll=res$stock[3],
				proportion=res$proportion,
				pi=res$p,
				mu=res$referenceRank,
				tik=res$tik,
				partition=res$cluster,
				entropy=res$entropy,
				probability=res$probability,
				convergence=TRUE,
				partial=FALSE,
				distanceZ=res$distZ,
				distanceMu=res$distMu,
				distanceProp=res$distProp,
				distancePi=res$distP,
				piInitial=res$initPi,
   				muInitial=res$initMu,
				proportionInitial=res$initProportion)
		}

		if(detail)
		{
			cat("RESULTS:\n")
			cat("NUMBER OF CLUSTERS: ",g)
			cat("\nLoglikelihood =",res$stock[3])   
			cat("\nBIC=",res$stock[4])
			cat("\nICL=",res$stock[5])		
			cat("\nProportion:",res$proportion)
			cat("\nProbabilities pi:\n")
			print(res$p)
			cat("\nReference ranks mu:\n")
			print(res$referenceRank)
		}
	}
	else
	{
		
		result=new(Class="Output",convergence=FALSE)
	}

	return(result)
}

