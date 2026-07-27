/*
 * math_utils.h
 *
 * Probability and statistical utility functions: LogSumExp, log-probability
 * normalization, Gibbs acceptance, BIC, and Rand index.
 */

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <vector>
#include <Eigen/Dense>


/**
 * Compute the Bayesian Information Criterion.
 * BIC = -2 * loglikelihood + nbParam * log(nbData)
 */
double BIC(double loglikelihood, int nbDonnees, int nbParam);

/**
 * Compute the Rand index between two partitions.
 * https://en.wikipedia.org/wiki/Rand_index
 */
double computeRandIndex(std::vector<int> const &z1, std::vector<int> const &z2);

/**
 * Metropolis-Hastings acceptance step for a Gibbs sampler.
 * Samples x ~ U(0,1) and accepts the candidate if x < p2 / (p1 + p2).
 *
 * @param logP1 log-probability of the current state
 * @param logP2 log-probability of the candidate
 * @return true if the candidate is accepted
 */
bool acceptChange(double logP1, double logP2);

/**
 * LogSumExp:
 *   LSE(x1, ..., xn) = log(exp(x1) + ... + exp(xn))
 *
 * Computed in a numerically stable way as:
 *   max(xi) + log( sum_i exp(xi - max(xi)) )
 */
double LSE(Eigen::ArrayXd &logProba);

/**
 * Normalize a vector of log-probabilities into plain probabilities.
 * Given (log p1, ..., log pn), returns (p1/S, ..., pn/S) where S = sum(pi).
 */
Eigen::ArrayXd normalizeLogProba(Eigen::ArrayXd &logProba);

/** In-place version of normalizeLogProba. */
void normalizeLogProbaInPlace(Eigen::ArrayXd &logProba);


#endif /* MATH_UTILS_H_ */
