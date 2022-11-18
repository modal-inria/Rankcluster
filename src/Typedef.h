#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <vector>
#include "Eigen/Dense"

typedef std::vector<int> Rank;
typedef std::vector<std::vector<int>> MultivariateRank;


struct PartialRank
{
  // observed rank
  Rank x;
  // order of presentation
  Rank y;
  // if true, the rank contains partial or ties
  bool isPartial;
  // missing element of the rank or ties
  std::vector<std::vector<int>> missingData;
  // index of the 0 or ties
  std::vector<std::vector<int>> missingIndex;
};

typedef std::vector<PartialRank> MultivariatePartialRank;

struct MuList
{
  int freq;
  std::vector<std::vector<Rank>> fullRank;
  std::vector<std::vector<double>> p;
  std::vector<double> prop;
  MuList *nextMu;
};


struct SEMparameters
{
    // number of iteration in the gibbs of SE step for each dimension of rank
    std::vector<int> nGibbsSE;
    // number of iteration in the gibbs of M step for each dimension of rank
    std::vector<int> nGibbsM;
    // maximum number of iteration of the SEM algorithm
    int maxIt;
    // burn-in period of SEM algorithm
    int burnAlgo;
    // number of iteration in the gibbs of the likelihood computation
    int nGibbsL;
    // burn-in period of the likelihood computation
    int burnL;
    // maximum number of try of the SEM
    int maxTry;
    // if true print information
    bool verbose;
};

struct OutParameters
{
  // log-likelihood
  double L;  //TODO rename into ll
  // bic criterion
  double bic;
  // icl criterion
  double icl;
  //
  Eigen::ArrayXXd tik; //TODO rename into posteriorProbabilities
  Eigen::ArrayXd entropy;
  //
  Eigen::ArrayXXd probabilities; //TODO rename into conditionalProbabilities
  // percentage of confidence in final estimation of missing data
  std::vector<std::vector<std::vector<double>>> partialRankScore;

  // algorithm initialization
  std::vector<std::vector<Rank>> initialPartialRank;
  std::vector<std::vector<double>> initialP;
  std::vector<int> initialZ;
  std::vector<double> initialProportion;
  std::vector<std::vector<Rank>> initialMu;

  // distance between parameters
  std::vector<std::vector<double>> distProp;
  std::vector<std::vector<std::vector<double>>> distP;
  std::vector<std::vector<std::vector<int>>> distMu;
  std::vector<double> distZ;
  std::vector<std::vector<std::vector<int>>> distPartialRank;
};

#endif // TYPEDEF_H