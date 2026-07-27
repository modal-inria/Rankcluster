/*
 * ISRfunctions.cpp
 *
 * Functions associated with the ISR (Insertion Sorting Rank) distribution.
 */

#include <algorithm>  // std::min_element
#include <cmath>      // std::log, std::exp
#include <cstddef>    // std::size_t
#include <map>
#include <Rmath.h>    // runif

#include "ISRfunctions.h"
#include "functions.h"

/*
 * Computation of A(x,y) (total comparisons) and G(x,y,mu) (good comparisons).
 *
 * Uses precomputed inverse-rank lookup tables to replace the repeated O(m)
 * linear scans in positionRank with O(1) array lookups.
 *
 * See appendix A and B: https://hal.archives-ouvertes.fr/hal-00441209v3/document
 */
ComparisonStats comparaison(Rank const &x, Rank const &y, Rank const &mu)
{
    int const m(static_cast<int>(mu.size()));

    // Build inverse-rank lookup tables (1-indexed elements, 0-based positions)
    std::vector<int> posX(m + 1, 0);
    for (int i = 0; i < m; i++)
        posX[x[i]] = i;

    std::vector<int> posMu(m + 1, 0);
    for (int i = 0; i < m; i++)
        posMu[mu[i]] = i;

    int gplus(0), gmoins(0);
    std::vector<int> ajmoins, ajplus;        // cumulative unions of Aj- and Aj+
    std::vector<int> ajplusb, ajmoinsb;      // per-j Aj+ and Aj- sets
    std::vector<int> ajplusbIndex;           // original indices for ajplusb
    ajplusb.reserve(m);
    ajmoinsb.reserve(m);
    ajplusbIndex.reserve(m);
    // ajplus has at most one entry per j iteration (m-1 total)
    ajplus.reserve(static_cast<std::size_t>(m - 1));
    // ajmoins can accumulate up to m*(m-1)/2 entries in the worst case
    ajmoins.reserve(static_cast<std::size_t>(m) * static_cast<std::size_t>(m - 1) / 2);

    for (int j = 1; j < m; j++)
    {
        for (int i = 0; i < j; i++)
        {
            if (posX[y[i]] < posX[y[j]])
            {
                // y[i] comes before y[j] in x  →  belongs to Aj-
                ajmoins.push_back(i);
                ajmoinsb.push_back(i);
            }
            else
            {
                // y[i] comes after y[j] in x  →  candidate for Aj+
                ajplusb.push_back(posX[y[i]]);
                ajplusbIndex.push_back(i);
            }
        }

        if (!ajplusb.empty())
        {
            // The representative of Aj+ is the element of {y[i]: i in ajplus_candidates}
            // that appears first in x, i.e. has the smallest position in x.
            int minIdx = static_cast<int>(
                std::min_element(ajplusb.begin(), ajplusb.end()) - ajplusb.begin());
            ajplus.push_back(ajplusbIndex[minIdx]);

            // Good comparison for Aj+: mu ranks y[j] before the Aj+ representative
            if (posMu[y[j]] < posMu[y[ajplus.back()]])
                gplus++;

            ajplusb.clear();
            ajplusbIndex.clear();
        }

        if (!ajmoinsb.empty())
        {
            // Good comparisons for Aj-: mu ranks each y[i] (i in Aj-) before y[j]
            for (int i : ajmoinsb)
            {
                if (posMu[y[i]] < posMu[y[j]])
                    gmoins++;
            }
            ajmoinsb.clear();
        }
    }

    return ComparisonStats{static_cast<int>(ajmoins.size() + ajplus.size()),
                           gmoins + gplus};
}

/** Compute log p(x|y; mu, p) under the ISR model */
double lnProbaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p)
{
    ComparisonStats cs = comparaison(x, y, mu);

    // Guard the two cases that produce NaN via 0*log(0); let all other cases
    // (including p=1 with total≠good, or p=0 with good≠0) produce -inf
    // naturally from the formula, matching the original behavior.
    if (p == 1.0 && cs.total == cs.good)
        return 0.;
    if (p == 0.0 && cs.good == 0)
        return 0.;

    return cs.good * std::log(p) + (cs.total - cs.good) * std::log(1. - p);
}

double probaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p)
{
    return std::exp(lnProbaCond(x, y, mu, p));
}


// Simulate n samples from an ISR(mu, p)
std::vector<Rank> simulISR(int const &n, int const &m, Rank const &mu, double const &p)
{
    std::vector<Rank> simul(n, Rank(m, 0));
    Rank s(m);
    initializeRank(s);

    // Build inverse-rank lookup for mu once (mu is constant across all draws)
    std::vector<int> posMu(m + 1, 0);
    for (int k = 0; k < m; k++)
        posMu[mu[k]] = k;

    for (int i = 0; i < n; i++)
    {
        // Draw a random presentation order s
        Rshuffle(s.begin(), s.end());
        simul[i][0] = s[0];

        for (int j = 1; j < m; j++)
        {
            int l = 0;
            bool inserted = false;
            while (!inserted && l < j)
            {
                double rnd = runif(0., 1.);
                bool goodComp = (posMu[s[j]] < posMu[simul[i][l]]);
                if ((goodComp && rnd < p) || (!goodComp && rnd > p))
                {
                    // Insert s[j] at position l, shifting elements right
                    for (int k = j - 1; k >= l; k--)
                        simul[i][k + 1] = simul[i][k];
                    simul[i][l] = s[j];
                    inserted = true;
                }
                else
                {
                    l++;
                }
            }
            if (!inserted)
                simul[i][j] = s[j];
        }
    }
    return simul;
}


// Simulate a sample from a mixture of ISR distributions
void simulMixtureISR(std::vector<Rank> &simul, std::vector<Rank> const &mu,
                     std::vector<double> const &p, std::vector<double> const &prop)
{
    int n = static_cast<int>(simul.size());
    int m = static_cast<int>(mu[0].size());

    Eigen::ArrayXd probaEig(prop.size());
    for (int k = 0; k < static_cast<int>(prop.size()); k++)
        probaEig(k) = prop[k];

    for (int i = 0; i < n; i++)
    {
        int classe = sampleMultinomial(probaEig);
        simul[i] = simulISR(1, m, mu[classe], p[classe])[0];
    }
}


double proba(std::vector<Rank> const &x, std::vector<Rank> const &mu,
             std::vector<double> const &pi)
{
    int d = static_cast<int>(pi.size());
    std::vector<double> probaDim(d, 0.);
    std::map<int, std::vector<int>> diffDim;

    for (int dim = 0; dim < d; dim++)
        diffDim[static_cast<int>(mu[dim].size())].push_back(dim);

    // Iterate over distinct rank sizes to avoid recomputing y for each dimension
    for (auto const &kv : diffDim)
    {
        int m = kv.first;
        std::vector<int> tabFact = tab_factorial(m);
        std::vector<int> listeY  = listIndexOrderOfPresentation(m, tabFact);
        double mult = 2. / static_cast<double>(tabFact[m - 1]);

        for (int i = 0; i < static_cast<int>(tabFact[m - 1]) / 2; i++)
        {
            Rank y = index2rankNoCheck(listeY[i], m, tabFact);
            for (int j : kv.second)
                probaDim[j] += probaCond(x[j], y, mu[j], pi[j]);
        }

        for (int j : kv.second)
            probaDim[j] *= mult;
    }

    double finalProba = probaDim[0];
    for (int dim = 1; dim < d; dim++)
        finalProba *= probaDim[dim];

    return finalProba;
}
