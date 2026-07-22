#include "RinterfaceTest.h"
#include "RinterfaceFunctions.h"

using namespace Rcpp;
using namespace std;

vector<vector<double>> convertToVVd(NumericMatrix data)
{
    int n = data.nrow(), p = data.ncol();
    vector<vector<double>> output(n, vector<double>(p));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            output[i][j] = data[j * n + i];

    return output;
}

vector<vector<int>> convertToVVi(NumericMatrix data)
{
    int n = data.nrow(), p = data.ncol();
    vector<vector<int>> output(n, vector<int>(p));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            output[i][j] = data[j * n + i];

    return output;
}

vector<RankStruct> downUniVariateRank(NumericMatrix XR)
{
    int const n(XR.nrow());
    int const m(XR.ncol());
    vector<RankStruct> donnees(n);
    set<int> element;

    for (int i(1); i < m + 1; i++)
        element.insert(i);

    for (int j(0); j < n; j++)
        donnees[j].rank.resize(m);

    for (int j(0); j < n; j++)
    {
        donnees[j].missingNumber = element;

        donnees[j].isPartial = false;

        for (int i(0); i < m; i++)
        {
            donnees[j].rank[i] = XR[j + i * n];
            if (donnees[j].rank[i] == 0)
            {
                donnees[j].isPartial = true;
                donnees[j].missingIndex.push_back(i);
            }
            else
                donnees[j].missingNumber.erase(donnees[j].rank[i]);
        }
    }

    return donnees;
}

// [[Rcpp::export(name = ".kullback")]]
double kullback(std::vector<int> m, NumericMatrix mu1, NumericMatrix mu2,
                NumericMatrix p1, NumericMatrix p2,
                std::vector<double> proportion1, std::vector<double> proportion2)
{
    vector<vector<double>> p1C = convertToVVd(p1);
    vector<vector<double>> p2C = convertToVVd(p2);
    vector<vector<vector<int>>> mu1C = numMat2vvvInt(mu1, m);
    vector<vector<vector<int>>> mu2C = numMat2vvvInt(mu2, m);

    return divKL(m, mu1C, mu2C, p1C, p2C, proportion1, proportion2);
}

// [[Rcpp::export]]
double adkhi2(NumericMatrix donnees, std::vector<double> p, std::vector<double> proportion,
              NumericMatrix mu, int nBootstrap)
{
    vector<vector<int>> muC = convertToVVi(mu);
    vector<vector<int>> data = convertToVVi(donnees);

    return khi2(data, p, proportion, muC, nBootstrap);
}

// [[Rcpp::export]]
double adkhi2partial(NumericMatrix donnees, std::vector<double> p, std::vector<double> proportion,
                     NumericMatrix mu, int nBootstrap)
{
    vector<vector<int>> muC = convertToVVi(mu);
    vector<RankStruct> data = downUniVariateRank(donnees);

    return khi2partial(data, p, proportion, muC, nBootstrap);
}

