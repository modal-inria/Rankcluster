/*
 * sampling_utils.cpp
 *
 * Random number generation, shuffling, and sampling utilities.
 */

#include <vector>
#include <Rmath.h>

#include "sampling_utils.h"

int randWrapper(const int n)
{
    return static_cast<int>(unif_rand() * n);
}

int sampleMultinomial(Eigen::ArrayXd const &proba)
{
    int g = static_cast<int>(proba.size());
    std::vector<double> lim(g + 1, 0.);

    for (int k = 0; k < g; k++)
        lim[k + 1] = lim[k] + proba(k);

    double alea = static_cast<double>(unif_rand());

    for (int j = 0; j < g; j++)
        if (lim[j] <= alea && alea < lim[j + 1])
            return j;

    return g - 1;
}
