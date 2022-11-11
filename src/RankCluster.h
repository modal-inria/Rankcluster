#ifndef RANKCLUSTER_H_
#define RANKCLUSTER_H_

/* @file RankCluster.h
 * @brief Implementation of methods of the class @c RankCluster
 * See https://hal.archives-ouvertes.fr/hal-00441209v3/document and
 * https://hal.inria.fr/hal-00743384/document for mathematical background
 */
//#include <RcppEigen.h>
#include <vector>
#include <set>
#include <utility>
#include "Eigen/Dense"

#include "ISRfunctions.h"
#include "functions.h"
#include "Typedef.h"


class RankCluster
{
public:
  // default constructor
  RankCluster();
  /**
   * @param X data one row= a multivariate rank
   * @param g number of clusters
   * @param m size of rank of each dimension
   * @param param parameters of SEM algorithm
   */
  RankCluster(std::vector<std::vector<int>> const &X, int g, std::vector<int> const &m, SEMparameters const &param);
  // constructor with given initialization of parameters
  RankCluster(std::vector<std::vector<int>> const &X, std::vector<int> const &m, SEMparameters const &param,
              std::vector<double> const &proportion, std::vector<std::vector<double>> const &p,
              std::vector<std::vector<std::vector<int>>> const &mu);

  // copy constructor
  RankCluster(RankCluster &rankClusterObject);

  // destructor
  virtual ~RankCluster();

  // run the SEM algorithm
  void run();

  // getters
  inline int d() const { return d_; }
  inline int n() const { return n_; }
  inline int g() const { return g_; }
  inline std::vector<int> m() const { return m_; }
  inline std::vector<int> z() const { return z_; }
  inline std::vector<std::vector<double>> p() const { return p_; }
  inline std::vector<std::vector<Rank>> mu() const { return mu_; }
  inline std::vector<double> proportion() const { return proportion_; }
  inline std::vector<std::vector<int>> indexPartialData() const { return indexPartialData_; }
  inline std::vector<std::vector<PartialRank>> data() const { return data_; }
  inline std::vector<int> rank(int dim, int index) const { return data_[dim][index].x; }
  inline bool dataOk() const { return dataOk_; }
  inline bool convergence() const { return convergence_; }
  inline bool partial() const { return partial_; }
  inline std::vector<std::vector<int>> indexPb() const { return indexPb_; }
  inline SEMparameters parameter() const { return parameter_; }

  // output getters
  inline Eigen::ArrayXXd tik() const { return output_.tik; }
  inline Eigen::ArrayXd entropy() const { return output_.entropy; }
  inline Eigen::ArrayXXd probabilities() const { return output_.probabilities; }
  Eigen::ArrayXd probability() const;
  inline double bic() const { return output_.bic; }
  inline double icl() const { return output_.icl; }
  inline double L() const { return output_.L; }
  inline std::vector<std::vector<Rank>> initialPartialRank() const { return output_.initialPartialRank; }
  inline std::vector<std::vector<double>> initialP() const { return output_.initialP; }
  inline std::vector<int> initialZ() const { return output_.initialZ; }
  inline std::vector<std::vector<Rank>> initialMu() const { return output_.initialMu; }
  inline std::vector<double> initialProportion() const { return output_.initialProportion; }
  inline std::vector<std::vector<double>> distProp() const { return output_.distProp; }
  inline std::vector<std::vector<std::vector<double>>> distP() const { return output_.distP; }
  inline std::vector<std::vector<std::vector<int>>> distMu() const { return output_.distMu; }
  inline std::vector<double> distZ() const { return output_.distZ; }
  inline std::vector<std::vector<std::vector<int>>> distPartialRank() const { return output_.distPartialRank; }
  inline std::vector<std::vector<std::vector<double>>> partialRankScore() const { return output_.partialRankScore; }

  /** re-estimation of criterion */
  void estimateCriterion(double &L, double &bic, double &icl);

protected:
  /** convert X in vector<vector<PartialRank>>
   * @param X raw data one row= a multi dimensional rank
   */
  void conversion2data(std::vector<std::vector<int>> const &X);

  /** read rank. used in conversion2data
   * @param X raw data one row= a multi dimensional rank
   * @param dim actual dimension
   * @param j actual index of the sample
   * @param indM transformation of m_
   */
  void readRankingRank(std::vector<std::vector<int>> const &X, int const &dim, int const &j, std::vector<int> const &indM);

  // initialization of parameters
  void initialization();
  void initializeZ();
  void initializeP();
  void initializeMu();
  void initializeY();
  void initializePartialRank();
  void fillIndexPartialData();
  void saveInitialization();

  void estimateProportion();


  /** SE step */
  void SEstep();

  /** univariate gibbs sampler for y estimation
   * @param indexDim index of the dimension (<d_)
   */
  void gibbsY(int indexDim);

  /** Compute tik (probability to belong to the different clusters) for the individual ind*/
  Eigen::ArrayXd computeTik(int const ind);

  /** simulate partition using the current estimated parameters */
  void sampleZ();

  /** univariate gibbs sampler for partial rank estimation
   * @param indexDim index of the dimension (<d_)
   */
  void gibbsX(int indexDim);

  /** M step */
  void Mstep();

  /** Estimate mu and p during the M step
   * @param indexDim index of the dimension
   * @param classNumber index of the cluster
   */
  void estimateMuP(int indexDim, int classNumber);

  /** Simulate a new candidate for mu
   * @param indexDim index of the dimension
   * @param classNumber index of the cluster
   * @param mu current mu (will be modified)
   * @param logPcurrent log-probability associated with mu (will be modified)
   */
  void simulateCandidateMuKJ(int indexDim, int classNumber, Rank &mu, double &logPcurrent);

  /** Update P parameter for the given dimension and class
   * @param indexDim index of the dimension
   * @param classNumber index of the cluster
   * @param mu current mu
   * @param sumG used in computeCompletedLoglikehoodKJ function (will be modified)
   * @param sumA_G used in computeCompletedLoglikehoodKJ function (will be modified)
   */
  double updatePKJ(int indexDim, int classNumber, int nObsClass, Rank const& mu, double &sumG, double &sumA_G);

  /** Compute the completed log-likelihood
   *
   * @param p current p parameter
   * @param sumG sum of good comparisons (computed in updatePKJ function)
   * @param sumA_G sum of all comparisons minus good comparisons (computed in updatePKJ function)
   */
  double computeCompletedLoglikehoodKJ(double p, double sumG, double sumA_G);


  /** Select the best set of parameters (the ones maximizing the log-likelihood)
   *
   * @param resMu mu for every iteration of the SEM
   * @param resP pi for every iteration of the SEM
   * @param resProp proportion for every iteration of the SEM
   */
  void selectBestParameters(std::vector<std::vector<std::vector<Rank>>> &resMu,
                            std::vector<std::vector<std::vector<double>>> &resP,
                            std::vector<std::vector<double>> &resProp);

  /** List the different found mu
   *
   * @param resMu mu for every iteration of the SEM
   * @param resP pi for every iteration of the SEM
   * @param resProp proportion for every iteration of the SEM
   */
  MuList * findDifferentMu(std::vector<std::vector<std::vector<Rank>>> &resMu,
                           std::vector<std::vector<std::vector<double>>> &resP,
                           std::vector<std::vector<double>> &resProp);

  /** Mean the parameters for a mu */
  void meanParameters(MuList * currMu);

  /** compute the log likelihood for a set of parameter
   * @param mu estimated central rank
   * @param p estimated p (dispersion parameter)
   * @param proportion estimated proportion
   * @param tik used for store tik
   * @param Y used for store estimated y
   * @param xTemp used for store estimated partial ranks
   * @param probabilities probability of each individual at each iteration
   * @param score used for confidence in estimated partial rank
   */
  double computeLogLikelihood(std::vector<std::vector<Rank>> const &mu,
                              std::vector<std::vector<double>> const &p,
                              std::vector<double> const &proportion, Eigen::ArrayXXd &tik,
                              std::vector<std::vector<Rank>> &Y,
                              std::vector<std::vector<Rank>> &xTemp, Eigen::ArrayXXd &probabilities,
                              std::vector<std::vector<std::vector<double>>> &score);

  /** compute the final partition (z_) */
  void computePartition();

  /** Store the parameters */
  void storeParameters(int const iterNumber, std::vector<std::vector<double>> &resProp,
                       std::vector<std::vector<std::vector<double>>> &resP,
                       std::vector<std::vector<std::vector<Rank>>> &resMu,
                       std::vector<std::vector<int>> &resZ,
                       std::vector<std::vector<std::vector<Rank>>> &resPartialData);

  /** compute distance between final parameters and each iteration parameters */
  void computeDistance(std::vector<std::vector<double>> const &resProp,
                       std::vector<std::vector<std::vector<double>>> const &resP,
                       std::vector<std::vector<std::vector<Rank>>> const &resMu,
                       std::vector<std::vector<int>> const &resZ,
                       std::vector<std::vector<std::vector<Rank>>> const &resDonneesPartiel);

private:
  // contains the size of rank for each dim
  std::vector<int> m_;
  // number of individuals
  int n_;
  // number of dimension
  int d_;
  // number of cluster
  int g_;
  // data of the form data[dimension][index]
  std::vector<std::vector<PartialRank>> data_;
  // estimated cluster of each individual
  std::vector<int> z_;
  // estimated rank parameter of each cluster:  mu_[dimension][cluster]
  std::vector<std::vector<Rank>> mu_;
  // estimated probability parameter of each cluster :  p_[dimension][cluster]
  std::vector<std::vector<double>> p_;
  // estimated proportion of the mixture model
  std::vector<double> proportion_;
  // algorithm parameters
  SEMparameters parameter_;
  // distance and initialization of the algorithm
  OutParameters output_;
  // true if there is partial rank in the data
  bool partial_;
  // index of partial data
  std::vector<std::vector<int>> indexPartialData_;
  // if true, SEM has converged
  bool convergence_;
  // if true, good data
  bool dataOk_;
  // index of rank with problem for each dimension
  std::vector<std::vector<int>> indexPb_;
};

#endif /* RANKCLUSTER_H_ */
