/*
 * sampling_utils.h
 *
 * Random number generation, shuffling, and sampling utilities.
 */

#ifndef SAMPLING_UTILS_H_
#define SAMPLING_UTILS_H_

#include <utility>      // std::swap
#include <Eigen/Dense>


/** Generate a random integer in [0, n-1] using R's RNG */
int randWrapper(const int n);

/**
 * Equivalent of std::shuffle backed by R's RNG.
 * Uses Fisher-Yates (Knuth) algorithm.
 */
template <class RandomAccessIterator>
void Rshuffle(RandomAccessIterator first, RandomAccessIterator last)
{
    for (auto i = (last - first) - 1; i > 0; --i)
        std::swap(first[i], first[randWrapper(static_cast<int>(i + 1))]);
}

/**
 * Sample one index according to a discrete distribution.
 * @param proba vector of non-negative probabilities (need not sum to 1)
 * @return sampled index in [0, proba.size()-1]
 */
int sampleMultinomial(Eigen::ArrayXd const &proba);


#endif /* SAMPLING_UTILS_H_ */
