/*
 * rank_utils.h
 *
 * Utilities for rank manipulation: position lookup, index/rank conversions,
 * Kendall distance, and related helpers.
 */

#ifndef RANK_UTILS_H_
#define RANK_UTILS_H_

#include <vector>
#include "Typedef.h"


/** Return the 0-based position of element i in rank x */
int positionRank(Rank const &x, int const &i);

/** Fill rank with 1, 2, ..., m */
void initializeRank(Rank &rank);

/** Fill rank with a random permutation of 1, 2, ..., m */
void randomRank(Rank &rank);

/**
 * Factorial (recursive).
 * @param nombre a positive integer
 * @return nombre!
 */
int factorial(int const &nombre);

/**
 * Compute all factorials from 1! to m!.
 * @param m a positive integer
 * @return vector of size m with tab[i] = (i+1)!
 */
std::vector<int> tab_factorial(int const &m);

/**
 * Convert a rank to its Lehmer-code index (1-based).
 * @param rank    a valid permutation rank
 * @param tabFactorial  precomputed factorials from tab_factorial()
 * @return index in [1, m!]
 */
int rank2index(Rank const &rank, std::vector<int> const &tabFactorial);

/** Convert a vector of ranks to their indices. */
std::vector<int> rank2index(std::vector<Rank> const &rankList,
                             std::vector<int> const &tabFact);

/**
 * Convert a 1-based Lehmer-code index to the corresponding rank.
 * Returns a zero-filled rank if index > m!.
 */
Rank index2rank(int index, int const &m, std::vector<int> const &tabFactorial);
/** Overload: computes tabFactorial internally. */
Rank index2rank(int index, int const &m);

/**
 * Like index2rank but skips the bounds check for performance.
 * The caller must guarantee index <= m!.
 */
Rank index2rankNoCheck(int index, int const &m,
                        std::vector<int> const &tabFactorial);

/**
 * Return the subset of order-of-presentation indices for which
 * computing p(x|y; mu, p) is necessary (exploiting the symmetry
 * p(x|y; mu, p) == p(x|swap(y[0],y[1]); mu, p)).
 *
 * @param m             size of the rank
 * @param tabFactorial  precomputed factorials from tab_factorial()
 * @return vector of length m!/2 containing 1-based indices
 */
std::vector<int> listIndexOrderOfPresentation(
    int const &m, std::vector<int> const &tabFactorial);

/**
 * Reverse a rank in-place (1 2 3 4 → 4 3 2 1).
 */
void invertRank(Rank &rank);

/**
 * Convert a rank from ordering representation to ranking representation.
 * The two representations are inverses of each other.
 */
Rank ordering2ranking(Rank const &x);

/**
 * Compute the Kendall tau distance between two ranks given in ordering
 * representation.
 * https://en.wikipedia.org/wiki/Kendall_tau_distance
 */
int distanceKendall(Rank const &x, Rank const &y);


#endif /* RANK_UTILS_H_ */
