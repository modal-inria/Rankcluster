/*
 * File containing functions associated with the ISR distribution
 */

#ifndef ISR_FUNCTIONS_H_
#define ISR_FUNCTIONS_H_

#include <vector>
#include "Typedef.h"


/**
 * compute  A(x,y) and G(x,y,mu)
 * @param x rank
 * @param y order of presentation of x
 * @param mu order of reference
 * @return a vector of 2 elements (A(x,y),G(x,y,mu))
 *
 * A(x,y)=total number of comparison in the insertion sorting algorithm
 * G(x,y,mu)= total number of good comparison according to mu in the insertion sorting algorithm
 *
 */
std::vector<int> comparaison(Rank const &x, Rank const &y, Rank const &mu);

/**
 * compute the conditional probability  p(x|y;mu,p)
 * @param x rank
 * @param y order of presentation of x
 * @param mu order of reference
 * @param p probability of making a good comparison
 * @return p(x|y;mu,p)
 */
double probaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p);
double lnProbaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p);


/**
 * compute probability of x according to multivariate ISR
 * @param x rank for each dimension for compute probability
 * @param mu reference rank for each dimension
 * @param pi dispersion parameter for each dimension
 * @return p(x;mu,pi)
 */
double proba(std::vector<Rank> const &x, std::vector<Rank> const &mu, std::vector<double> const &pi);


/**
 * simulation of a n-sample of ISR(mu,p)
 * @param n size of the sample
 * @param m size of the rank
 * @param mu rank
 * @param p probability of a good comparison
 * @return a n-sample of ISR(mu,p)
 */
std::vector<Rank> simulISR(int const &n, int const &m, Rank const &mu, double const &p);

/**
 * Simulate a sample of mixture of ISR
 * @param simul sample will be modified
 * @param mu reference rank
 * @param p dispersion parameter: probability of a good comparison
 * @param prop proportion of the mixture
 */
void simulMixtureISR(std::vector<Rank> &simul, std::vector<Rank> const &mu, std::vector<double> const &p,
                     std::vector<double> const &prop);


#endif /* ISR_FUNCTIONS_H_ */