#include "RinterfaceSEM.h"

using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
List semR(NumericMatrix X, std::vector<int> m, int K, int Qsem, int Bsem,
          int Ql, int Bl, std::vector<int> RjSE, std::vector<int> RjM,
          int maxTry, int run, bool verbose)
{
    SEMparameters param;
    param.nGibbsSE = RjSE;
    param.nGibbsM = RjM;
    param.maxIt = Qsem;
    param.burnAlgo = Bsem;
    param.nGibbsL = Ql;
    param.burnL = Bl;
    param.maxTry = maxTry;
    param.verbose = verbose;

    int n(X.nrow()), col(X.ncol());
    vector<vector<int>> data(n, vector<int>(col));
    for (int i(0); i < n; i++)
    {
        for (int j(0); j < col; j++)
        {
            data[i][j] = X[i + j * n];
        }
    }

    RankCluster semgibbs(data, K, m, param);
    if (!semgibbs.dataOk())
    {
        vector<double> stock(1, 2);
        return List::create(Named("stock") = wrap(stock), Named("indexPb") = wrap(semgibbs.indexPb()));
    }

    semgibbs.run();

    //multiple run
    if (run > 1)
    {
        RankCluster semgibbsb(data, K, m, param);
        double L(-numeric_limits<double>::max());
        if (semgibbs.convergence())
            L = semgibbs.L();

        for (int i(1); i < run; i++)
        {
            semgibbsb.run();
            if (semgibbsb.convergence())
            {
                if (semgibbsb.L() > L)
                {
                    L = semgibbsb.L();
                    semgibbs = semgibbsb;
                }
            }
        }
    }

    if (semgibbs.convergence())
    {
        vector<double> stock(11);
        stock[0] = 1;
        stock[1] = semgibbs.partial();
        stock[2] = semgibbs.L();
        stock[3] = semgibbs.bic();
        stock[4] = semgibbs.icl();
        int d = m.size();
        vector<vector<vector<int>>> rankData(d, vector<vector<int>>(n));
        vector<vector<vector<double>>> scorePart(d, vector<vector<double>>(n));

        for (int i(0); i < n; i++)
        {
            for (int dim(0); dim < d; dim++)
            {
                rankData[dim][i] = semgibbs.rank(dim, i);
                scorePart[dim][i].resize(m[dim], 1);
            }
        }

        vector<vector<vector<int>>> dataInit(rankData);

        if (semgibbs.partial())
        {
            vector<vector<vector<int>>> dataInitSEM(semgibbs.initialPartialRank());
            vector<vector<int>> indexPartial = semgibbs.indexPartialData();
            vector<vector<vector<double>>> scorePartial(semgibbs.partialRankScore());

            for (int dim(0); dim < d; dim++)
            {
                int compteur(0);
                for (vector<int>::iterator it = indexPartial[dim].begin(); it != indexPartial[dim].end(); it++)
                {
                    dataInit[dim][*it] = dataInitSEM[dim][compteur];
                    scorePart[dim][*it] = scorePartial[dim][compteur];
                    compteur++;
                }
            }
        }

        return List::create(
            //parameters
            Named("stock") = wrap(stock),
            Named("referenceRank") = wrap(semgibbs.mu()),
            Named("p") = wrap(semgibbs.p()),
            Named("proportion") = wrap(semgibbs.proportion()),
            Named("cluster") = wrap(semgibbs.z()),
            Named("tik") = wrap(semgibbs.tik()),
            Named("entropy") = wrap(semgibbs.entropy()),
            Named("partialRank") = wrap(rankData),
            Named("probability") = wrap(semgibbs.probability()),
            //distance
            Named("distP") = wrap(semgibbs.distP()),
            Named("distMu") = wrap(semgibbs.distMu()),
            Named("distProp") = wrap(semgibbs.distProp()),
            Named("distZ") = wrap(semgibbs.distZ()),
            Named("distPartialRank") = wrap(semgibbs.distPartialRank()),
            //initialization
            Named("initMu") = wrap(semgibbs.initialMu()),
            Named("initZ") = wrap(semgibbs.initialZ()),
            Named("initPi") = wrap(semgibbs.initialP()),
            Named("initPartialRank") = wrap(dataInit),
            Named("initProportion") = wrap(semgibbs.initialProportion()),
            Named("scorePartial") = wrap(scorePart)
        );
    }
    else
    {
        vector<double> stock(1, 0);
        return List::create(Named("stock") = wrap(stock));
    }
}
