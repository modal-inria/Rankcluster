#include "functions.h"
#include "runFunctions.h"
#include "runTest.h"
#include "RankCluster.h"

using namespace Rcpp ;
using namespace std ;
using namespace Eigen ;

vector<vector<vector<int> > > numMat2vvvInt(NumericMatrix XR,vector<int> const& m)
{
    int indiceElement(0);
    int const g(XR.nrow()),d(m.size());
    vector<vector<vector<int> > > donnees(d,vector<vector<int> > (g));

    //taille des lignes du fichier à importer
    vector<int> indM(d+1,0);
    for(int i(0);i<d;i++)
        indM[i+1]=indM[i]+m[i];


    for(int i(0);i<d;i++)
        for(int j(0);j<g;j++)
            donnees[i][j].resize(m[i]);

    for(int j(0);j<g;j++)
    {
        for(int k(0);k<d;k++)
        {
            indiceElement=0;

            for(int i(indM[k]);i<indM[k+1];i++)
            {
                donnees[k][j][indiceElement]=XR[j+i*g];
                indiceElement++;
            }

        }

    }

    return donnees;
}


RcppExport SEXP freqMultiR(SEXP X,SEXP m)
{
	NumericVector mR(m);
	vector<int> M=as<vector<int> > (mR);
	NumericMatrix XR(X);
	int const n(XR.cols()),d(M.size());
	vector<vector<vector<int> > > donnees(d,vector<vector<int> > (n));

	donnees=numMat2vvvInt(XR,M);

	pair<vector<vector<vector<int> > >,vector<int> > res;
	res=freqMulti(donnees);

	int taille(0),compteur(0);
	for(int i(0);i<d;i++)
		taille+=res.first[i][0].size();

	vector<vector<int> > data(res.first[0].size(),vector<int>(taille));

	for(int j(0);j<data.size();j++)
	{
		compteur=0;
		for(int i(0);i<d;i++)
		{
			for(int k(0);k<M[i];k++)
			{
				data[j][compteur]=res.first[i][j][k];
				compteur++;
			}
		}
	}

	return List::create(Named("data")=wrap(data),Named("freq")=wrap(res.second));

}

RcppExport SEXP simulISRR(SEXP n,SEXP m,SEXP mu,SEXP p)
{
	NumericVector muR(mu);
	vector<int> muC=as<vector<int> > (muR);
	int const nC=as<int>(n),mC=as<int>(m);
	double const pC=as<double>(p);

	vector<vector<int> > simul;
	simul=simulISR(nC,mC,muC,pC);

	NumericMatrix data(nC,mC);
	for(int i(0);i<nC;i++)
		for(int j(0);j<mC;j++)
			data(i,j)=simul[i][j];


	return data;
}

RcppExport SEXP loglikelihood(SEXP X,SEXP mu,SEXP p, SEXP proportion,SEXP m, SEXP iterL, SEXP burnL)
{
    //conversion
	NumericMatrix XR(X);
	int n(XR.nrow()),col(XR.ncol());
	vector<vector<int> > data(n,vector<int> (col));
	for(int i(0);i<n;i++)
		for(int j(0);j<col;j++)
			data[i][j]=XR[i+j*n];

    NumericVector proportionR(proportion);
	NumericVector mR(m);
    vector<int> mC=as<vector<int> > (mR);
	vector<double> prop=as<vector<double> > (proportionR);
	vector<vector<double> > pC;
	pC=convertToVVd(p);

	NumericMatrix muR(mu);
	vector<vector<vector<int> > > muC;
	muC=numMat2vvvInt(mu,mC);

    SEMparameters param;
	param.nGibbsSE=mC;
	param.nGibbsM=mC;
	param.maxIt=1;
	param.burnAlgo=1;
	param.nGibbsL=as<int>(iterL);
	param.burnL=as<int>(burnL);
	param.maxTry=1;
	param.detail=false;

    RankCluster estimLog(data,mC,param,prop,pC,muC);

    double L,bic,icl;
    estimLog.estimateCriterion(L,bic,icl);

    return List::create(Named("ll")=wrap(L),Named("bic")=wrap(bic),Named("icl")=wrap(icl));
}
