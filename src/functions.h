/*
 * functions.h
 *
 *  Created on: 1 mars 2013
 *      Author: Quentin Grimonprez
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <Rmath.h>

#include "Typedef.h"


/**r random number generator, code from : http://gallery.rcpp.org/articles/r-function-from-c++/
 * generate an integer between 0 and n - 1
 */
int randWrapper(const int n);

/** equivalent of std::shuffle with an R generator */
template <class RandomAccessIterator>
void Rshuffle(RandomAccessIterator first, RandomAccessIterator last)
{
  for (auto i = (last - first) - 1; i > 0; --i)
  {
    std::swap(first[i], first[randWrapper(i + 1)]);
  }
}

/** Initialize a rank as 1 2 ... m */
void initializeRank(Rank &rank);

/** Initialize a rank as a random rank */
void randomRank(Rank &rank);

/**
 * search the position of i in x
 * @param x rank
 * @param i integer which we want the position in x
 * @return the position of i in x
 */
int positionRank(Rank const &x, int const &i);

/**
 * factorial function (recursive)
 * @param m an integer >0
 * @return m!
 */
int factorial(int const &nombre);

/**
 * compute all the factorial from 1! to m!
 * @param m a positive integer
 * @return a vector of size m containing all the factorial from 1! to m!
 */
std::vector<int> tab_factorial(int const &m);

/**
 *-----------rank->index------------
 * convert a rank into an integer
 * @param rang rank
 * @param  tabFactorial \see tab_factorial
 * @return index
 */
int rank2index(Rank const &rank, std::vector<int> const &tabFactorial);

/**
 *-----------rank->index------------
 * @see rank2index for a vector of rank
 * @param rankList vector of rank
 * @param  tabFactorial @see tab_factorial
 * @return vector of index
 */
std::vector<int> rank2index(std::vector<Rank> const &rankList, std::vector<int> const &tabFact);

/**
 * conversion index->rank
 * @param index index of the rank
 * @param m size of the rank
 * @param tabFactorial vector containing 1! to m! (@see tab_factorial)(optional)
 * @return the rank corresponding to the index
 * index | rank
 * 0 | 1 2 ... m
 *
 * m! | m m-1 ... 1
 */
Rank index2rank(int index, int const &m, std::vector<int> const &tabFactorial);
Rank index2rankNoCheck(int index, int const &m, std::vector<int> const &tabFactorial);
Rank index2rank(int index, int const &m);

/**
 * Return a vector containing the index of order of presentation y for which we have to compute probability
 * @param m  size of the rank
 * @param tabFactorial vector containing 1! to m! (@see tab_factorial)
 * @return a vector containing the index of order of presentation y for which we have to compute probability
 *
 * p(x|y;mu,p) doesn't change if the 2 first element of y are inverted
 * So we compute probabilities for the half of y
 */
std::vector<int> listIndexOrderOfPresentation(int const &m, std::vector<int> const &tabFactorial);

/**
 * invert a rank (1 2 3 become 3 2 1)
 * @param rank rank to invert
 */
void invertRank(Rank &rank);

/**
 * Compute the BIC
 * @param loglikelihood the loglikelihood of the data
 * @param nbDonnees total number of data
 * @return BIC
 */
double BIC(double loglikelihood, int nbDonnees, int nbParam);

/**
 * Compute the Rand index between 2 partitions
 * @param z1 partition
 * @param z2 partition
 * @return a double, the Rand index
 */
double computeRandIndex(std::vector<int> const &z1, std::vector<int> const &z2);

/**
 * Conversion from ordering representation to ranking representation
 *
 * @param x a rank
 * @return a rank: the ranking representation of x
 */
Rank ordering2ranking(Rank const &x);

/**
 *  Compute the Kendall's distance between 2 ranks (in ordering representation)
 *
 * https://en.wikipedia.org/wiki/Kendall_tau_distance
 *
 * @param x a rank in ordering representation
 * @param y a rank in ordering representation of the same size than x
 * @return an integer, the kendall distance between x and y
 */
int distanceKendall(Rank const &x, Rank const &y);

/**
 * Sort the parameters such that the first cluster is the cluster with the more little index of mu
 * @param mu index of the rank of the first dimension of listeMu
 * @param p parameter of the ISR
 * @param prop proportion of the mixture model
 * @param listeMu reference rank
 * @param z partition
 * @param g number of cluster
 * @param d number of dimension
 * @param n number of individual
 *
 * listeMu, p, prop and z are modified if necessary
 *
 */
void tri_insertionMulti(Rank &mu, std::vector<double> &prop, std::vector<std::vector<double>> &p,
                        std::vector<std::vector<Rank>> &listeMu, std::vector<int> &z, int const &g,
                        int const &d, int const &n);

/**
 * compute frequency of a data set
 * @param rankList data
 * @return a pair with the unique rank and the frequency of each unique rank
 */
std::pair<std::vector<std::vector<Rank>>, std::vector<int>> freqMulti(std::vector<std::vector<Rank>> const &rankList);


/**
 * accept or not a change in a gibbs sampler
 *
 * Sample x~U(0, 1), accept the change if x < p2/(p1 + p2)
 *
 * @param logP1 log-probability of the current element
 * @param logP2 log-probability of the candidate
 * @return true if the change is accepted
 */
bool acceptChange(double const logP1, double const logP2);

/** LogSumExp function
 * https://en.wikipedia.org/wiki/LogSumExp
 *
 * LSE(x1, x2, ..., xn) = (max xi) + log(exp(x1 - (max xi)) + ... + exp(xn - (max xi)))
 * LSE(x1, x2, ..., xn) = log(exp(x1) + ... + exp(xn))
 *
 * Application for proba: LSE(log(p1), ..., log(pn)) = log(p1 + ... + pn)
 */
double LSE(Eigen::ArrayXd &logProba);

/**
 * LSE trick to normalize a vector of log probabilities
 * Given (log(p1), ..., log(pn)), it returns (p1/sum(pi), ..., pn/sum(pi))
 *
 * Useful to compute tik
 */
Eigen::ArrayXd normalizeLogProba(Eigen::ArrayXd &logProba);
void normalizeLogProbaInPlace(Eigen::ArrayXd &logProba);

/**
 * Sample according to a multinomial with given probabilities
 *
 * @param proba vector of probabilities of each element
 */
int sampleMultinomial(Eigen::ArrayXd const &proba);

#endif /* FUNCTIONS_H_ */
