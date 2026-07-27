/*
 * functions.cpp
 *
 *  Created on: 1 mars 2013
 *      Author: Quentin Grimonprez
 *
 * Higher-level cluster utilities: parameter sorting and frequency counting.
 * Low-level primitives have been moved to rank_utils, math_utils, and
 * sampling_utils.
 */

#include <map>

#include "functions.h"

// mu: index of first-dimension elements per cluster
// listeMu: listMu[dim][cluster][elem]
void tri_insertionMulti(Rank &mu, std::vector<double> &prop,
                        std::vector<std::vector<double>> &p,
                        std::vector<std::vector<Rank>> &listeMu,
                        std::vector<int> &z,
                        int const &g, int const &d, int const &n)
{
    int elem, elemorder;
    double elemprop;
    std::vector<double> elemp(d);
    std::vector<std::vector<int>> elemmu(d);
    std::vector<int> order(g);
    for (int l = 0; l < g; l++)
        order[l] = l;

    // Insertion sort on mu[0..g-1]; replicate swaps across all dimensions
    for (int i = 1; i < g; i++)
    {
        elem      = mu[i];
        elemprop  = prop[i];
        elemorder = order[i];
        for (int l = 0; l < d; l++)
        {
            elemp[l]  = p[l][i];
            elemmu[l] = listeMu[l][i];
        }

        int j = i;
        for (; j > 0 && mu[j - 1] > elem; j--)
        {
            order[j]  = order[j - 1];
            mu[j]     = mu[j - 1];
            prop[j]   = prop[j - 1];
            for (int l = 0; l < d; l++)
            {
                p[l][j]      = p[l][j - 1];
                listeMu[l][j] = listeMu[l][j - 1];
            }
        }
        order[j]  = elemorder;
        mu[j]     = elem;
        prop[j]   = elemprop;
        for (int l = 0; l < d; l++)
        {
            p[l][j]      = elemp[l];
            listeMu[l][j] = elemmu[l];
        }
    }

    // Re-order the partition z according to the cluster reordering
    for (int l = 0; l < n; l++)
    {
        for (int k = 0; k < g; k++)
        {
            if (z[l] == order[k])
            {
                z[l] = k;
                break;
            }
        }
    }
}


/*
 * Compute the frequency table of a multivariate rank dataset.
 *
 * Replaces the original hand-rolled linked-list implementation with a
 * std::map lookup, preserving insertion order of first-seen ranks.
 */
std::pair<std::vector<std::vector<Rank>>, std::vector<int>>
freqMulti(std::vector<std::vector<Rank>> const &rankList)
{
    const int d(static_cast<int>(rankList.size()));
    const int N(static_cast<int>(rankList[0].size()));

    // Map from multivariate rank key → index into donnees/freq
    // (preserves insertion order of first occurrence)
    std::map<std::vector<Rank>, int> indexMap;
    std::vector<std::vector<Rank>> donnees(d);
    std::vector<int> freq;

    for (int j = 0; j < N; j++)
    {
        std::vector<Rank> key(d);
        for (int i = 0; i < d; i++)
            key[i] = rankList[i][j];

        auto it = indexMap.find(key);
        if (it != indexMap.end())
        {
            freq[it->second]++;
        }
        else
        {
            int idx = static_cast<int>(freq.size());
            indexMap[key] = idx;
            for (int i = 0; i < d; i++)
                donnees[i].push_back(key[i]);
            freq.push_back(1);
        }
    }

    return std::make_pair(donnees, freq);
}
