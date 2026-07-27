/*
 * ISRfunctions.h
 *
 * Functions associated with the ISR (Insertion Sorting Rank) distribution.
 * See https://hal.archives-ouvertes.fr/hal-00441209v3/document for the
 * mathematical background.
 */

#ifndef ISR_FUNCTIONS_H_
#define ISR_FUNCTIONS_H_

#include <vector>
#include "Typedef.h"


/**
 * Result of the comparison-counting step of the ISR algorithm.
 *
 * A(x,y)       = total number of comparisons performed by the insertion
 *                sorting algorithm when sorting x according to presentation
 *                order y.
 * G(x,y,mu)    = number of those comparisons that agree with the reference
 *                rank mu.
 */
struct ComparisonStats
{
    int total;  ///< A(x,y)
    int good;   ///< G(x,y,mu)
};


/**
 * Compute A(x,y) and G(x,y,mu).
 * See appendix A and B of https://hal.archives-ouvertes.fr/hal-00441209v3
 *
 * @param x   rank (observed ordering)
 * @param y   order of presentation of x
 * @param mu  reference rank
 * @return ComparisonStats with total = A(x,y) and good = G(x,y,mu)
 */
ComparisonStats comparaison(Rank const &x, Rank const &y, Rank const &mu);

/**
 * Compute log p(x|y; mu, p) under the ISR model.
 * @param x   rank
 * @param y   order of presentation of x
 * @param mu  reference rank
 * @param p   probability of making a good comparison
 * @return log p(x|y; mu, p)
 */
double lnProbaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p);

/**
 * Compute p(x|y; mu, p) = exp(lnProbaCond(x, y, mu, p)).
 */
double probaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p);


/**
 * Compute the marginal probability p(x; mu, pi) of a multivariate rank x
 * under the multivariate ISR model.
 *
 * @param x   multivariate rank (one rank per dimension)
 * @param mu  reference ranks (one per dimension)
 * @param pi  dispersion parameters (one per dimension)
 * @return p(x; mu, pi)
 */
double proba(std::vector<Rank> const &x, std::vector<Rank> const &mu,
             std::vector<double> const &pi);


/**
 * Simulate n independent draws from ISR(mu, p).
 * @param n  number of samples
 * @param m  size of the rank
 * @param mu reference rank
 * @param p  probability of a good comparison
 * @return n-by-m matrix of sampled ranks
 */
std::vector<Rank> simulISR(int const &n, int const &m, Rank const &mu,
                            double const &p);

/**
 * Simulate a sample from a mixture of ISR distributions.
 * @param simul  output sample (modified in-place, size determines n)
 * @param mu     reference ranks, one per component
 * @param p      dispersion parameters, one per component
 * @param prop   mixture proportions
 */
void simulMixtureISR(std::vector<Rank> &simul, std::vector<Rank> const &mu,
                     std::vector<double> const &p,
                     std::vector<double> const &prop);


#endif /* ISR_FUNCTIONS_H_ */