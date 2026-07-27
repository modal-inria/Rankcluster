/*
 * functions.h
 *
 * Backward-compatible facade: includes all sub-modules and declares the
 * remaining higher-level functions (tri_insertionMulti, freqMulti).
 *
 * New code should prefer including the focused sub-headers directly:
 *   - sampling_utils.h  (randWrapper, Rshuffle, sampleMultinomial)
 *   - rank_utils.h      (positionRank, rank2index, index2rank, …)
 *   - math_utils.h      (LSE, normalizeLogProba, acceptChange, BIC, …)
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "sampling_utils.h"
#include "rank_utils.h"
#include "math_utils.h"

#include <vector>
#include <utility>
#include "Typedef.h"


/**
 * Sort cluster parameters so that the first cluster has the smallest
 * first-dimension mu index (insertion sort on mu[0]).
 *
 * @param mu     Lehmer indices of the first dimension's mu, one per cluster
 * @param prop   mixture proportions (modified in-place)
 * @param p      dispersion parameters p[dim][cluster] (modified in-place)
 * @param listeMu  reference ranks listeMu[dim][cluster] (modified in-place)
 * @param z      partition vector z[individual] (modified in-place)
 * @param g      number of clusters
 * @param d      number of dimensions
 * @param n      number of individuals
 */
void tri_insertionMulti(Rank &mu, std::vector<double> &prop,
                        std::vector<std::vector<double>> &p,
                        std::vector<std::vector<Rank>> &listeMu,
                        std::vector<int> &z,
                        int const &g, int const &d, int const &n);

/**
 * Compute the frequency table of a multivariate rank dataset.
 *
 * @param rankList  data as rankList[dimension][observation]
 * @return pair of (unique_ranks[dimension][unique_index], counts)
 *         in order of first occurrence
 */
std::pair<std::vector<std::vector<Rank>>, std::vector<int>>
freqMulti(std::vector<std::vector<Rank>> const &rankList);


#endif /* FUNCTIONS_H_ */
