/*
 * math_utils.cpp
 *
 * Probability and statistical utility functions.
 */

#include <cmath>
#include <Rmath.h>

#include "math_utils.h"

double BIC(double loglikelihood, int nbDonnees, int nbParam)
{
    return -2.0 * loglikelihood + nbParam * std::log(static_cast<double>(nbDonnees));
}

double computeRandIndex(std::vector<int> const &z1, std::vector<int> const &z2)
{
    const int N(static_cast<int>(z1.size()));
    double a(0.), b(0.), c(0.), d(0.);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            bool sameZ1 = (z1[i] == z1[j]);
            bool sameZ2 = (z2[i] == z2[j]);
            if (sameZ1 && sameZ2)       a++;
            else if (!sameZ1 && !sameZ2) b++;
            else if (sameZ1 && !sameZ2)  c++;
            else                         d++;
        }
    }

    return (a + b) / (a + b + c + d);
}

bool acceptChange(double const logP1, double const logP2)
{
    double logSum;
    if (logP1 > logP2)
        logSum = logP1 + std::log(std::exp(logP2 - logP1) + 1.0);
    else
        logSum = logP2 + std::log(std::exp(logP1 - logP2) + 1.0);

    double proba2 = std::exp(logP2 - logSum);
    return (runif(0., 1.) < proba2);
}

double LSE(Eigen::ArrayXd &logProba)
{
    double m = logProba.maxCoeff();
    return m + std::log((logProba - m).exp().sum());
}

Eigen::ArrayXd normalizeLogProba(Eigen::ArrayXd &logProba)
{
    double lse = LSE(logProba);
    return (logProba - lse).exp();
}

void normalizeLogProbaInPlace(Eigen::ArrayXd &logProba)
{
    double lse = LSE(logProba);
    logProba = (logProba - lse).exp();
}
