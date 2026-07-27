#include "functions.h"
#include "RinterfaceFunctions.h"
#include "RinterfaceTest.h"
#include "RankCluster.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace std;
using namespace Eigen;

vector<vector<vector<int>>> numMat2vvvInt(NumericMatrix XR, vector<int> const &m)
{
    int indiceElement(0);
    int const g(XR.nrow()), d(m.size());
    vector<vector<vector<int>>> donnees(d, vector<vector<int>>(g));

    //taille des lignes du fichier à importer
    vector<int> indM(d + 1, 0);
    for (int i(0); i < d; i++)
        indM[i + 1] = indM[i] + m[i];

    for (int i(0); i < d; i++)
        for (int j(0); j < g; j++)
            donnees[i][j].resize(m[i]);

    for (int j(0); j < g; j++)
    {
        for (int k(0); k < d; k++)
        {
            indiceElement = 0;

            for (int i(indM[k]); i < indM[k + 1]; i++)
            {
                donnees[k][j][indiceElement] = XR[j + i * g];
                indiceElement++;
            }
        }
    }

    return donnees;
}

// [[Rcpp::export]]
List freqMultiR(NumericMatrix X, std::vector<int> m)
{
    int const n(X.nrow()), d(m.size());
    vector<vector<vector<int>>> donnees(d, vector<vector<int>>(n));

    donnees = numMat2vvvInt(X, m);

    pair<vector<vector<vector<int>>>, vector<int>> res;
    res = freqMulti(donnees);

    int taille(0), compteur(0);
    for (int i(0); i < d; i++)
        taille += res.first[i][0].size();

    vector<vector<int>> data(res.first[0].size(), vector<int>(taille));

    for (int j(0); j < (int)data.size(); j++)
    {
        compteur = 0;
        for (int i(0); i < d; i++)
        {
            for (int k(0); k < m[i]; k++)
            {
                data[j][compteur] = res.first[i][j][k];
                compteur++;
            }
        }
    }

    return List::create(Named("data") = wrap(data), Named("freq") = wrap(res.second));
}

// [[Rcpp::export]]
NumericMatrix simulISRR(int n, int m, std::vector<int> mu, double p)
{
    vector<vector<int>> simul = simulISR(n, m, mu, p);

    NumericMatrix data(n, m);
    for (int i(0); i < n; i++)
        for (int j(0); j < m; j++)
            data(i, j) = simul[i][j];

    return data;
}

// [[Rcpp::export]]
List loglikelihood(NumericMatrix X, NumericMatrix mu, NumericMatrix p,
                   std::vector<double> proportion, std::vector<int> m,
                   int iterL, int burnL, int IC, int nb_cpus)
{
    int n(X.nrow()), col(X.ncol());
    vector<vector<int>> data(n, vector<int>(col));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < col; j++)
            data[i][j] = X[i + j * n];

    vector<vector<double>> pC = convertToVVd(p);
    vector<vector<vector<int>>> muC = numMat2vvvInt(mu, m);

    SEMparameters param;
    param.nGibbsSE = m;
    param.nGibbsM = m;
    param.maxIt = 1;
    param.burnAlgo = 1;
    param.nGibbsL = iterL;
    param.burnL = burnL;
    param.maxTry = 1;
    param.verbose = false;

    RankCluster estimLog(data, m, param, proportion, pC, muC);
    if (!estimLog.dataOk())
    {
        return List::create(Named("ll") = wrap("pb"));
    }

    vector<double> L(IC, 0), bic(IC, 0), icl(IC, 0);

#ifdef _OPENMP
    nb_cpus = std::min(nb_cpus, omp_get_num_procs());
#pragma omp parallel num_threads(nb_cpus)
    {
#pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < IC; i++)
        {
            RankCluster estimLogpar(estimLog);
            estimLogpar.estimateCriterion(L[i], bic[i], icl[i]);
        }
    }
#else
    {
        for (int i = 0; i < IC; i++)
        {
            RankCluster estimLogpar(estimLog);
            estimLogpar.estimateCriterion(L[i], bic[i], icl[i]);
        }
    }
#endif

    return List::create(Named("ll") = wrap(L),
                        Named("bic") = wrap(bic),
                        Named("icl") = wrap(icl));
}

// [[Rcpp::export]]
double computeProba(NumericMatrix X, NumericMatrix mu, std::vector<double> pi, std::vector<int> m)
{
    vector<vector<int>> muC(m.size());
    vector<vector<int>> xC(m.size());

    for (int i = 0; i < static_cast<int>(m.size()); i++)
    {
        xC[i].resize(m[i]);
        muC[i].resize(m[i]);
    }

    int j = 0, k = 0;
    for (int i = 0; i < mu.ncol(); i++)
    {
        if (k == m[j])
        {
            j++;
            k = 0;
        }
        muC[j][k] = mu[i];
        xC[j][k] = X[i];

        k++;
    }

    return proba(xC, muC, pi);
}
