/* @file RankCluster.cpp
 * @brief Implementation of methods of the class @c RankCluster
 * See https://hal.archives-ouvertes.fr/hal-00441209v3/document and
 * https://hal.inria.fr/hal-00743384/document for mathematical background
 */

#include <cstdlib>
#include <algorithm>
#include <string>
#include <limits>
#include <cmath>
#include <map>
#include <iostream>
#include <ctime>
#include <Rmath.h>

#include "RankCluster.h"

using namespace std;
using namespace Eigen;
//using namespace Rcpp;


/**********************************************************************************************
 *
 * CONSTRUCTORS, DESTRUCTOR
 *
 **********************************************************************************************/

RankCluster::RankCluster()
{
}

RankCluster::RankCluster(std::vector<std::vector<int>> const &X, int g,
                         vector<int> const &m, SEMparameters const &param)
  : m_(m), n_(X.size()), d_(m.size()), g_(g),
    data_(d_, vector<PartialRank>(n_)),
    z_(n_),
    mu_(d_, vector<Rank>(g_)),
    p_(d_, vector<double>(g_)),
    proportion_(g),
    parameter_(param),
    partial_(false),
    dataOk_(true),
    indexPb_(m.size())
{
  //  try
  //  {
  // convert data in the good notation and create information about missing and partial
  conversion2data(X);
  //  }
  //  catch(string const& errorString)
  //    {dataOk_=false;}
}

RankCluster::RankCluster(vector<vector<int>> const &X, vector<int> const &m, SEMparameters const &param,
                         vector<double> const &proportion, vector<vector<double>> const &p,
                         vector<vector<Rank>> const &mu)
  : m_(m), n_(X.size()), d_(m.size()), g_(proportion.size()),
    data_(d_, vector<PartialRank>(n_)),
    z_(n_), mu_(mu), p_(p),
    proportion_(proportion),
    parameter_(param),
    partial_(false),
    dataOk_(true),
    indexPb_(m.size())
{
  //  try
  //  {
  // convert data in the good notation and create information about missing and partial
  conversion2data(X);
  //  }
  //  catch(string const& errorString)
  //    {dataOk_=false;}
}

// copy constructor
RankCluster::RankCluster(RankCluster &rankClusterObject)
  : m_(rankClusterObject.m()),
    n_(rankClusterObject.n()),
    d_(rankClusterObject.d()),
    g_(rankClusterObject.g()),
    data_(rankClusterObject.data()),
    mu_(rankClusterObject.mu()),
    p_(rankClusterObject.p()),
    proportion_(rankClusterObject.proportion()),
    parameter_(rankClusterObject.parameter()),
    partial_(rankClusterObject.partial()),
    dataOk_(rankClusterObject.dataOk())
{
  // nothing to do
}

// destructor
RankCluster::~RankCluster()
{
  // nothing
}


/**********************************************************************************************
 *
 * READ DATA
 *
 **********************************************************************************************/

void RankCluster::readRankingRank(vector<vector<int>> const &X, int const &dim, int const &j, vector<int> const &indM)
{
  //initialization
  int indiceElement = 0;
  data_[dim][j].isPartial = false;

  //multi dim rank temporary
  vector<vector<int>> temp(m_[dim] + 1);

  for (int i = indM[dim]; i < indM[dim + 1]; i++)
  {
    temp[X[j][i]].push_back(indiceElement + 1);
    indiceElement++;
  }

  //vector containing index of partial element
  vector<int> partialIndex;

  int skip = 0;
  //index 0 is for missing, we don't manage in this loop
  for (int i = 1; i < (int)temp.size(); i++)
  {
    if (skip)
    {
      if (temp[i].size() != 0)
      {
        dataOk_ = false;
        indexPb_[dim].push_back(j + 1);
        //throw string("Problem with data.");
      }

      skip--;
    }
    else
    {
      //tied case
      if (temp[i].size() > 1)
      {
        data_[dim][j].isPartial = true;
        partial_ = true;
        skip = temp[i].size() - 1;
        data_[dim][j].missingData.push_back(temp[i]);

        vector<int> missingIndex(temp[i].size());

        for (int ii = 0; ii < (int)temp[i].size(); ii++)
          missingIndex[ii] = i + ii - 1;

        data_[dim][j].missingIndex.push_back(missingIndex);
      }
      else
      {
        //normal case
        if (temp[i].size() == 1)
          data_[dim][j].x[i - 1] = temp[i][0];
        else //temp[i].size=0//partial case
          partialIndex.push_back(i - 1);
      }
    }
  }

  //problem with the data : index of 0 et element at missing position don't match
  if (temp[0].size() != partialIndex.size())
  {
    dataOk_ = false;
    indexPb_[dim].push_back(j + 1);
    //throw string("Problem with data.");
  }

  //add partial
  if (temp[0].size() != 0)
  {
    data_[dim][j].isPartial = true;
    partial_ = true;
    data_[dim][j].missingData.push_back(temp[0]);
    data_[dim][j].missingIndex.push_back(partialIndex);
  }
}

void RankCluster::conversion2data(vector<vector<int>> const &X)
{
  //size of a row of X
  vector<int> indM(d_ + 1, 0);
  for (int i = 0; i < d_; i++)
    indM[i + 1] = indM[i] + m_[i];

  //resize data
  for (int i = 0; i < d_; i++)
    for (int j = 0; j < n_; j++)
      data_[i][j].x.resize(m_[i]);

  //begin the read of the data row by row
  for (int j = 0; j < n_; j++)
  {
    //dim by dim
    for (int dim(0); dim < d_; dim++)
    {
      //read rank j of dim dim
      readRankingRank(X, dim, j, indM);
    }
  }
}


/**********************************************************************************************
 *
 * INITIALIZATION METHODS
 *
 **********************************************************************************************/
void RankCluster::initialization()
{
  // t0 = clock();
  initializeZ();
  initializeP();
  initializeMu();
  estimateProportion();
  initializePartialRank();
  fillIndexPartialData();
  saveInitialization();
  // t1 = clock();

  // if(parameter_.verbose)
  //     cout << "Initialization completed in " <<(double) (t1-t0)/CLOCKS_PER_SEC << "s." << endl;
}

// random initialization of z_: multinomial law of size g_
void RankCluster::initializeZ()
{
  for (int i = 0; i < n_; i++)
    z_[i] = randWrapper(g_);
}

// initialization of p_ with double between 0.5 and 1
void RankCluster::initializeP()
{
  for (int k = 0; k < d_; k++)
  {
    for (int i = 0; i < g_; i++)
    {
      p_[k][i] = (double) runif(0.5, 1.);
    }
  }
}

// random initialization of mu_ with random rank of size m_
void RankCluster::initializeMu()
{
  for (int k = 0; k < d_; k++)
  {
    for (int i = 0; i < g_; i++)
    {
      mu_[k][i].resize(m_[k]);
      randomRank(mu_[k][i]);
    }
  }
}

// estimate proportion using z_: proportion of each class among z_
void RankCluster::estimateProportion()
{
  for (int k = 0; k < g_; k++)
    proportion_[k] = 0;

  for (int i = 0; i < n_; i++)
    proportion_[z_[i]]++;

  for (int k = 0; k < g_; k++)
    proportion_[k] /= (double)n_;
}

// order of presentation random initialization
void RankCluster::initializeY()
{
  for (int dim = 0; dim < d_; dim++)
  {
    vector<int> rankTemp(m_[dim]);
    initializeRank(rankTemp);
    for (int ind = 0; ind < n_; ind++)
    {
      Rshuffle(rankTemp.begin(), rankTemp.end());
      data_[dim][ind].y = rankTemp;
    }
  }
}

// order of partial rank: random initialization with missing index
void RankCluster::initializePartialRank()
{
  for (int dim = 0; dim < d_; dim++)
  {
    for (int ind = 0; ind < n_; ind++)
    {
      if (data_[dim][ind].isPartial)
      {
        for (int ii = 0; ii < (int)data_[dim][ind].missingIndex.size(); ii++)
        {
          Rank rankTemp(data_[dim][ind].missingIndex[ii]);
          Rshuffle(rankTemp.begin(), rankTemp.end());

          for (int iii = 0; iii < (int)data_[dim][ind].missingData[ii].size(); iii++)
            data_[dim][ind].x[rankTemp[iii]] = data_[dim][ind].missingData[ii][iii];
        }
      }
    }
  }
}


void RankCluster::fillIndexPartialData()
{
  indexPartialData_ = vector<vector<int>>(d_);
  for (int dim = 0; dim < d_; dim++)
  {
    for (int ind = 0; ind < n_; ind++)
    {
      if (data_[dim][ind].isPartial)
        indexPartialData_[dim].push_back(ind);
    }
  }
}

// save initialization in output_ member
void RankCluster::saveInitialization()
{
  vector<vector<Rank>> partialRankData(d_);
  for (int dim = 0; dim < d_; dim++)
    for (vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
      partialRankData[dim].push_back(data_[dim][*it].x);

  output_.initialPartialRank = partialRankData;
  output_.initialP = p_;
  output_.initialZ = z_;
  output_.initialMu = mu_;
  output_.initialProportion = proportion_;
}



/**********************************************************************************************
 *
 * SE STEP METHODS
 *
 **********************************************************************************************/

void RankCluster::SEstep()
{
  // simulation of order of presentation for each dimension
  for (int dim = 0; dim < d_; dim++)
    gibbsY(dim);

  sampleZ();

  // simulation of partial rank for each dimension
  for (int dim = 0; dim < d_; dim++)
    gibbsX(dim);
}

void RankCluster::gibbsY(int indexDim)
{
  double logPcurrent(0), logPcandidate(0);
  set<int>::iterator itset;

  //rank 1 2..m
  // Rank yTemp(m_[indexDim]);
  // initializeRank(yTemp);

  for (int ind = 0; ind < n_; ind++)
  {
    //Gibbs sampling
    Rank y(m_[indexDim]), yCandidate(m_[indexDim]), yCurrent(m_[indexDim]);

    //initialization of p1 and y1
    randomRank(y);
    yCurrent = y;
    logPcurrent = lnProbaCond(data_[indexDim][ind].x, yCurrent, mu_[indexDim][z_[ind]], p_[indexDim][z_[ind]]);

    //start of iteration
    for (int iter(0); iter < parameter_.nGibbsSE[indexDim]; iter++)
    {
      for (int i(0); i < m_[indexDim] - 1; i++)
      {
        // new y to test (old y with inversion of 2 adjacent elements)
        yCandidate = y;
        yCandidate[i] = y[i + 1];
        yCandidate[i + 1] = y[i];

        // compute the probability of accept the change of y
        logPcandidate = lnProbaCond(data_[indexDim][ind].x, yCandidate,
                                    mu_[indexDim][z_[ind]], p_[indexDim][z_[ind]]);

        bool changeAccepted = acceptChange(logPcurrent, logPcandidate);
        if (changeAccepted) //change acceptation
        {
          y = yCandidate;
          logPcurrent = logPcandidate;
          yCurrent = y;
        }
        else
          y = yCurrent;
      }
    }
    data_[indexDim][ind].y = y;
  }
}


Eigen::ArrayXd RankCluster::computeTik(int const ind)
{
  Eigen::ArrayXd tik = ArrayXd::Zero(g_);

  for (int k(0); k < g_; k++)
  {
    for (int l(0); l < d_; l++)
      tik(k) += lnProbaCond(data_[l][ind].x, data_[l][ind].y, mu_[l][k], p_[l][k]);

    tik[k] += log(proportion_[k]);
  }
  normalizeLogProbaInPlace(tik);

  return tik;
}

/* z follow a multinomial law of parameter tik */
void RankCluster::sampleZ()
{
  if (g_ != 1)
  {
    Eigen::ArrayXd tik(g_);

    for (int ind(0); ind < n_; ind++)
    {
      // computation of the probability to belong to each cluster
      tik = computeTik(ind);

      z_[ind] = sampleMultinomial(tik);
    }
  }
  else
  {
    for (int ind(0); ind < n_; ind++)
      z_[ind] = 0;
  }
}

void RankCluster::gibbsX(int indexDim)
{
  double logPcurrent(0), logPcandidate(0);

  for (int ind = 0; ind < n_; ind++)
  {
    if (data_[indexDim][ind].isPartial)
    {
      // Gibbs algorithm
      Rank x(m_[indexDim]), xCurrent(m_[indexDim]), xCandidate(m_[indexDim]);

      // mu and log-probability initialization
      x = data_[indexDim][ind].x;
      xCurrent = x;
      logPcurrent = lnProbaCond(xCurrent, data_[indexDim][ind].y, mu_[indexDim][z_[ind]], p_[indexDim][z_[ind]]);

      for (int iter = 0; iter < parameter_.nGibbsSE[indexDim]; iter++)
      {
        for (int ii = 0; ii < (int)data_[indexDim][ind].missingIndex.size(); ii++)
        {
          for (int i = 0; i < (int)data_[indexDim][ind].missingIndex[ii].size() - 1; i++)
          {
            // new rank to test: old rank where 2 partial elements are exchanged
            xCandidate = x;
            xCandidate[data_[indexDim][ind].missingIndex[ii][i]] = x[data_[indexDim][ind].missingIndex[ii][i + 1]];
            xCandidate[data_[indexDim][ind].missingIndex[ii][i + 1]] = x[data_[indexDim][ind].missingIndex[ii][i]];

            logPcandidate = lnProbaCond(xCandidate, data_[indexDim][ind].y, mu_[indexDim][z_[ind]],
                                        p_[indexDim][z_[ind]]);

            bool changeAccepted = acceptChange(logPcurrent, logPcandidate);
            if (changeAccepted)
            {
              x = xCandidate;
              logPcurrent = logPcandidate;
              xCurrent = x;
            }
            else
              x = xCurrent;
          }
        }
      }
      data_[indexDim][ind].x = x;
    }
  }
}


/**********************************************************************************************
 *
 * M STEP METHODS
 *
 **********************************************************************************************/

void RankCluster::Mstep()
{
  // update proportion
  estimateProportion();

  for (int k(0); k < g_; k++)
  {
    if (proportion_[k] == 0)
      throw string("Algorithm did not converge: a proportion is equal to 0");
  }

  // simulation of mu for each dim and cluster
  for (int dim(0); dim < d_; dim++)
  {
    for (int classNumber(0); classNumber < g_; classNumber++)
      estimateMuP(dim, classNumber);
  }
}



void RankCluster::estimateMuP(int indexDim, int classNumber)
{
  vector<int> tabFact = tab_factorial(m_[indexDim]);

  // vector to store results of each iteration
  vector<Rank> MU(parameter_.nGibbsM[indexDim]);
  vector<double> P(parameter_.nGibbsM[indexDim], 0.), L(parameter_.nGibbsM[indexDim], 0.);

  // elements will be grouped per mu. If the new mu is the same, we do not recompute P and likelihood
  map<int, vector<double>> resPerMu;
  map<int, vector<double>>::iterator iteratorResPerMu;
  int indMu;

  // initialization of mu
  Rank mu = mu_[indexDim][classNumber];

  // initial completed log-likelihood
  double logPcurrent = 0.;
  int sizeCluster = 0.;
  for (int ind(0); ind < n_; ind++)
  {
    if (z_[ind] == classNumber)
    {
      logPcurrent += lnProbaCond(data_[indexDim][ind].x, data_[indexDim][ind].y, mu, p_[indexDim][classNumber]);
      sizeCluster++;
    }
  }

  // Gibbs algorithm
  for (int iter(0); iter < parameter_.nGibbsM[indexDim]; iter++)
  {
    simulateCandidateMuKJ(indexDim, classNumber, mu, logPcurrent);
    MU[iter] = mu;

    indMu = rank2index(MU[iter], tabFact);
    iteratorResPerMu = resPerMu.find(indMu);

    if (iteratorResPerMu == resPerMu.end()) //if we have already tested this mu, we do not redo computation
    {
      double sumG(0.), sumA_G(0.);
      vector<double> vecTemp(2);

      P[iter] = updatePKJ(indexDim, classNumber, sizeCluster, mu, sumG, sumA_G);

      L[iter] = computeCompletedLoglikehoodKJ(P[iter], sumG, sumA_G);

      vecTemp[0] = P[iter];
      vecTemp[1] = L[iter];
      resPerMu[indMu] = vecTemp;
    }
    else
    {
      P[iter] = (iteratorResPerMu->second)[0];
      L[iter] = (iteratorResPerMu->second)[1];
    }

    p_[indexDim][classNumber] = P[iter];
  }

  // find the best (mu, p) according to completed log-likelihood
  int indice(max_element(L.begin(), L.end()) - L.begin());
  p_[indexDim][classNumber] = P[indice];
  mu_[indexDim][classNumber] = MU[indice];
}


void RankCluster::simulateCandidateMuKJ(int indexDim, int classNumber, Rank &mu, double &logPcurrent)
{
  Rank muCandidate, muCurrent(mu);  // ? is muCurrent necessary ?
  double logPcandidate;

  // new mu
  for (int k(0); k < m_[indexDim] - 1; k++)
  {
    // new mu to test
    muCandidate = mu;
    muCandidate[k] = mu[k + 1];
    muCandidate[k + 1] = mu[k];

    logPcandidate = 0.;

    // new probability
    for (int ind(0); ind < n_; ind++)
    {
      if (z_[ind] == classNumber)
      {
        logPcandidate += lnProbaCond(data_[indexDim][ind].x, data_[indexDim][ind].y,
                                     muCandidate, p_[indexDim][classNumber]);
      }
    }

    // acceptation of change or not
    bool changeAccepted = acceptChange(logPcurrent, logPcandidate);
    if  (changeAccepted)
    {
      mu = muCandidate;
      logPcurrent = logPcandidate;
      muCurrent = mu;
    }
    else
      mu = muCurrent;
  }

}

/*
 * sumG and sumA_G are modified to be used in computeCompletedLoglikehoodKJ function
 */
double RankCluster::updatePKJ(int indexDim, int classNumber, int nObsClass, Rank const& mu,
                              double &sumG, double &sumA_G)
{
  double s1 = 0.;
  vector<int> comp(2);

  // computation of G and A-G (number of good and bad comparison) for log-likelihood
  sumG = 0.;
  sumA_G = 0.;

  for (int ind(0); ind < n_; ind++)
  {
    if (z_[ind] == classNumber)
    {
      comp = comparaison(data_[indexDim][ind].x, data_[indexDim][ind].y, mu);

      s1 += comp[0];
      sumG += comp[1];
      sumA_G += (comp[0] - comp[1]);
    }
  }

  return sumG / s1;
}

/*
 * Compute the completed log-likelihood to choose the best couple (mu_k^j, p_k^j).
 *
 * Maximizing the completed log-likelihood for a new couple (mu_k^j, p_k^j)
 * (dimension j = indexDim, cluster k = classNumber)
 * is equivalent to maximizing
 * \sum_{individuals i of class k} P(x_i^k | y_i^k; mu_k^j, p_k^j) =
 * \sum_{individuals i of class k} G_i^j log(p_k^j) + \sum_{individuals i of class k} A_i^j - G_i^j) log(1 - p_k^j)
 */
double RankCluster::computeCompletedLoglikehoodKJ(double p, double sumG, double sumA_G)
{
  double loglikelihood = 0.;
  if ((p != 0) && (p != 1))
  {
    loglikelihood = sumG * log(p) + sumA_G * log(1 - p);
  }
  else
  {
    if ((p == 0) && (sumG == 0))
      loglikelihood = 0.;
    else
    {
      if ((p == 1) && (sumA_G == 0))
        loglikelihood = 0.;
      else
        loglikelihood = -numeric_limits<double>::max();
    }
  }

  return loglikelihood;
}

/**********************************************************************************************
 *
 * LIKELIHOOD, CRITERIA AND OTHER OUTPUTS METHODS
 *
 **********************************************************************************************/
void RankCluster::selectBestParameters(vector<vector<vector<Rank>>> &resMu,
                                       vector<vector<vector<double>>> &resP,
                                       vector<vector<double>> &resProp)
{
  // double t1, t2, tL(0);

  // t1 = clock();

  MuList *headMu = findDifferentMu(resMu, resP, resProp);
  MuList *currMu = headMu;
  MuList *next = 0;

  // t2 = clock();
  // cout << "Temps regroupement mu: " << (double) (t2-t1)/CLOCKS_PER_SEC << "s" << endl;

  // compute log-likelihood
  // if(parameter_.verbose)
  //     cout << "Number of reference rank which must compute the log-likelihood: " << nbMu << endl;

  double Llast(-numeric_limits<double>::max()), L;

  vector<vector<Rank>> Y(d_, vector<Rank>(n_)), xPartialTemp(output_.initialPartialRank);

  vector<vector<vector<double>>> scoreTemp(output_.initialPartialRank.size());
  for (int ii = 0; ii < (int)scoreTemp.size(); ii++)
  {
    scoreTemp[ii].resize(output_.initialPartialRank[ii].size());
    for (int iii = 0; iii < (int)scoreTemp[ii].size(); iii++)
    {
      scoreTemp[ii][iii].resize(output_.initialPartialRank[ii][iii].size());
    }
  }

  // Now, we have the list of all the different Mu
  currMu = headMu;
  ArrayXXd tik(n_, g_);
  ArrayXXd probabilities(n_, g_); // estimate probability for an individual to belong to each cluster
  bool hasNextElement;

  // for each mu, we will compute the associated log likelihood
  do
  {
    // if(parameter_.verbose)
    //     cout<<"*";

    meanParameters(currMu);

    // compute the log likelihood
    // t1 = clock();
    L = computeLogLikelihood(currMu->fullRank, currMu->p, currMu->prop, tik, Y, xPartialTemp, probabilities, scoreTemp);
    // t2 = clock();
    // tL += t2-t1;

    if (L > Llast)
    {
      // the current mu has a better log-likelihood, we save the parameter
      Llast = L;
      mu_ = currMu->fullRank;
      p_ = currMu->p;
      proportion_ = currMu->prop;
      output_.tik = tik;
      output_.L = L;
      output_.probabilities = probabilities;
      output_.partialRankScore = scoreTemp;
      for (int dim = 0; dim < d_; dim++)
      {
        for (int ind = 0; ind < n_; ind++)
          data_[dim][ind].y = Y[dim][ind];

        int compteur(0);
        for (vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
        {
          data_[dim][*it].x = xPartialTemp[dim][compteur];
          compteur++;
        }
      }
    }

    next = currMu->nextMu;
    hasNextElement = (currMu->nextMu != 0);
    delete currMu; // delete the mu
    currMu = next;
  }
  while(hasNextElement);


  // if(parameter_.verbose)
  // {
  //     cout << "Computing time for log-likelihood approximation: " << (double)tL / CLOCKS_PER_SEC
  //          << "s (" << (double)tL / CLOCKS_PER_SEC / compteur << "s per mu)." << endl;
  // }

}


MuList * RankCluster::findDifferentMu(vector<vector<vector<Rank>>> &resMu, vector<vector<vector<double>>> &resP,
                                      vector<vector<double>> &resProp)
{
  bool notFound(true), equalRank;

  MuList *headMu = new MuList;
  MuList *currMu = headMu;
  currMu->freq = 1;
  currMu->nextMu = 0;
  currMu->fullRank = resMu[0];
  currMu->p = resP[0];
  currMu->prop = resProp[0];

  int nbMu(1);

  for (int it(1); it < (parameter_.maxIt - parameter_.burnAlgo); it++) //we see all the mu
  {
    notFound = true;
    currMu = headMu;
    while (notFound)
    {
      equalRank = true;
      // look if the j-th mu is the same that the current mu
      for (int j(0); j < d_; j++)
      {
        for (int k(0); k < g_; k++)
        {
          for (int i(0); i < m_[j]; i++)
            if (currMu->fullRank[j][k][i] != resMu[it][j][k][i])
            {
              equalRank = false;
              break;
            }
        }
      }

      if (equalRank)
      {
        currMu->freq++;
        // we sum the proportion and p
        for (int k(0); k < g_; k++)
        {
          currMu->prop[k] += resProp[it][k];
          for (int j(0); j < d_; j++)
            currMu->p[j][k] += resP[it][j][k];
        }
        notFound = false; // no need to check other mu of the struct
      }
      else
      {
        // if the current mu is the last, we add the j-th mu in the struct
        if (currMu->nextMu == 0)
        {
          nbMu++;
          notFound = false;
          currMu->nextMu = new MuList;
          currMu = currMu->nextMu;
          currMu->freq = 1;
          currMu->nextMu = 0;
          currMu->fullRank = resMu[it];
          currMu->prop = resProp[it];
          currMu->p = resP[it];
        }
        else
          currMu = currMu->nextMu; // we will test the next mu
      }
    }
  }

  return headMu;
}


void RankCluster::meanParameters(MuList * currMu)
{
  for (int k(0); k < g_; k++)
  {
    currMu->prop[k] /= (double) currMu->freq;
    for (int j(0); j < d_; j++)
      currMu->p[j][k] /= (double) currMu->freq;
  }
}


// LL gibbs
double RankCluster::computeLogLikelihood(vector<vector<Rank>> const &mu, vector<vector<double>> const &p,
                                         vector<double> const &proportion, ArrayXXd &tik,
                                         vector<vector<Rank>> &Y, vector<vector<Rank>> &xTemp,
                                         ArrayXXd &probabilities,
                                         vector<vector<vector<double>>> &score)
{
    long double p1(0), p2(0), p1x(0), p2x(0), alea(0), l(0), li(0);
    double div((double)(parameter_.nGibbsL - parameter_.burnL));
    vector<int> compteur(d_, 0);
    Rank x1, x2;

    //objet pour stockage de calcul pour éviter répétition
    ArrayXXd proba1(d_, g_), proba2(d_, g_);
    ArrayXd proba1X(g_), proba2X(g_);

    //store proportion
    ArrayXd prop(g_);
    for (int i = 0; i < g_; i++)
        prop(i) = proportion[i];
    ArrayXXd propb(1, g_);
    propb = prop.transpose();

    vector<Rank> yTemp(d_), x(d_);
    vector<vector<vector<int>>> scoreCount(d_);
    for (int j = 0; j < d_; j++)
    {
        //génération rang 1 2 3 ..m par dimension
        yTemp[j].resize(m_[j]);
        for (int i = 0; i < m_[j]; i++)
            yTemp[j][i] = i + 1;

        //initialize score
        scoreCount[j].resize(m_[j], vector<int>(m_[j], 0));
    }

    //start estimation
    for (int ind = 0; ind < n_; ind++)
    {
        //cout<<"ind "<<ind<<endl;
        vector<Rank> y(d_), y2(d_), y1(d_);
        li = 0;
        //Gibbs Algorithm to sample yi

        //initialisation de y et p pour Gibbs
        y = yTemp;
        for (int j = 0; j < d_; j++)
        {
            for (int jj = 0; jj < m_[j]; jj++)
                for (int jjj = 0; jjj < m_[j]; jjj++)
                    scoreCount[j][jj][jjj] = 0;

            Rshuffle(y[j].begin(), y[j].end()); //permutation de 1 2 3 ..m
            x[j] = data_[j][ind].x;
        }

        y1 = y;

        for (int k = 0; k < g_; k++)
        {
            tik(ind, k) = 0;
            for (int j = 0; j < d_; j++)
                proba1(j, k) = probaCond(x[j], y1[j], mu[j][k], p[j][k]);
        }

        p1 = (long double)(propb * proba1.colwise().prod()).sum();
        proba2 = proba1;

        //start gibbs for sample ind
        for (int iter = 0; iter < parameter_.nGibbsL; iter++)
        {

            /*simulation des y*/
            for (int J = 0; J < d_; J++)
            {
                for (int K = 0; K < m_[J] - 1; K++)
                {
                    //"état" à tester (inversion de 2 éléments adjacents)
                    y2 = y;
                    y2[J][K] = y[J][K + 1];
                    y2[J][K + 1] = y[J][K];

                    //tester un stockage des proba calculées pour éviter répétition de calculs dans la boucle
                    for (int k = 0; k < g_; k++)
                        proba2(J, k) = probaCond(x[J], y2[J], mu[J][k], p[J][k]);

                    p2 = (long double)(propb * proba2.colwise().prod()).sum();

                    alea = (long double)runif(0., p1 + p2); //unif(0,p1+p2)(double)

                    if (alea < p2) //acceptation du changement de y
                    {
                        y[J] = y2[J];
                        p1 = p2;
                        proba1.row(J) = proba2.row(J);
                        y1[J] = y[J];
                    }
                    else //on ne modifie pas y
                    {
                        y[J] = y1[J]; //rajout J
                        proba2.row(J) = proba1.row(J);
                    }
                }
            }
            //y_i est mis à jour

            /*simulation des x_i^j qui sont partiels*/
            for (int J = 0; J < d_; J++)
            {
                if (data_[J][ind].isPartial) //simulation de xi si partiel
                {
                    x1 = x[J];
                    proba1X = proba1.row(J);
                    p1x = (proba1X * prop).sum();
                    for (int kk = 0; kk < (int)(data_[J][ind].missingIndex).size() - 1; kk++)
                    {
                        for (int k = 0; k < (int)(data_[J][ind].missingIndex[kk]).size() - 1; k++) //Gibbs sur les x
                        {
                            //nouveau x à tester
                            x2 = x[J];
                            x2[data_[J][ind].missingIndex[kk][k]] = x[J][data_[J][ind].missingIndex[kk][k + 1]];
                            x2[data_[J][ind].missingIndex[kk][k + 1]] = x[J][data_[J][ind].missingIndex[kk][k]];

                            for (int l = 0; l < g_; l++)
                                proba2X(l) = probaCond(x2, y[J], mu[J][l], p[J][l]);

                            p2x = (proba2X * prop).sum();

                            alea = (double)runif(0., p1x + p2x);

                            if (alea < p2x) //acceptation du changement
                            {
                                x[J] = x2;
                                p1x = p2x;
                                proba1X = proba2X;
                                x1 = x[J];
                            }
                            else
                                x[J] = x1;
                        }
                    }

                    proba1.row(J) = proba1X;
                }
            }

            if (iter >= parameter_.burnL)
            {
                ArrayXd calculInter(g_);
                for (int cl = 0; cl < g_; cl++)
                {
                    calculInter(cl) = 1;
                    for (int dim = 0; dim < d_; dim++)
                        calculInter(cl) *= proba1(dim, cl);

                    probabilities(ind, cl) = calculInter(cl);
                    calculInter(cl) *= propb(cl);
                }

                double den(calculInter.sum());
                tik.row(ind) += (calculInter / den);

                li += (long double)1 / den;

                //compute score partial rank
                for (int dim = 0; dim < d_; dim++)
                {
                    if (data_[dim][ind].isPartial)
                    {
                        for (int indElem = 0; indElem < m_[dim]; indElem++)
                            scoreCount[dim][indElem][x[dim][indElem] - 1]++;
                    }
                }
            } //end if not burn

        } //fin du gibbs pour l'individu ind

        //l -= log(li*div);
        l -= log(li);

        tik.row(ind) /= div;
        probabilities.row(ind) /= div;

        //sauvegarde des nouveau y et x
        for (int j(0); j < d_; j++)
        {
            Y[j][ind] = y[j];
            if (data_[j][ind].isPartial)
            {
                xTemp[j][compteur[j]] = x[j];
                for (int elem = 0; elem < m_[j]; elem++)
                    score[j][compteur[j]][elem] = ((double)(scoreCount[j][elem][x[j][elem] - 1]) / (double)(parameter_.nGibbsL - parameter_.burnL));
                compteur[j]++;
            }
        }

    } //fin boucle sur n

    l += n_ * log(div);
    return l;
}

/** compute the final partition
 * The final partition is the argmax per row of tik
 */
void RankCluster::computePartition()
{
  // TODO can be done directly with Eigen
  // see https://stackoverflow.com/questions/11430588/find-rowwise-maxcoeff-and-index-of-maxcoeff-in-eigen

  if (g_ > 1)
  {
    double max;
    for (int ind(0); ind < n_; ind++)
    {
      max = output_.tik(ind, 0);
      z_[ind] = 0;
      for (int k(1); k < g_; k++)
      {
        if (output_.tik(ind, k) > max)
        {
          max = output_.tik(ind, k);
          z_[ind] = k;
        }
      }
    }
  }
}


void RankCluster::computeDistance(vector<vector<double>> const &resProp, vector<vector<vector<double>>> const &resP,
                                  vector<vector<vector<Rank>>> const &resMu, vector<vector<int>> const &resZ,
                                  vector<vector<vector<Rank>>> const &resDonneesPartiel)
{
  int const iterTotal(parameter_.maxIt - parameter_.burnAlgo);

  //initialization of container
  output_.distProp = vector<vector<double>>(iterTotal, vector<double>(g_));
  output_.distP = vector<vector<vector<double>>>(iterTotal, vector<vector<double>>(d_, vector<double>(g_)));
  output_.distMu = vector<vector<vector<int>>>(iterTotal, vector<vector<int>>(d_, vector<int>(g_)));
  output_.distZ = vector<double>(iterTotal);

  //compute the distance between the output parameters and parameters from each iteration
  for (int i(0); i < iterTotal; i++)
  {
    //distance between partition
    output_.distZ[i] = computeRandIndex(z_, resZ[i]);

    for (int cl(0); cl < g_; cl++)
    {
      //distance between proportion
      output_.distProp[i][cl] = pow(resProp[i][cl] - proportion_[cl], 2);

      for (int dim(0); dim < d_; dim++)
      {
        //distance between p
        output_.distP[i][dim][cl] = pow(resP[i][dim][cl] - p_[dim][cl], 2);

        //distance between mu
        output_.distMu[i][dim][cl] = distanceKendall(mu_[dim][cl], resMu[i][dim][cl]);
      }
    }
  }

  //distance between partial rank
  vector<vector<vector<int>>> distRangPartiel(iterTotal, vector<vector<int>>(d_));
  if (partial_)
  {
    for (int i(0); i < iterTotal; i++)
    {
      for (int dim(0); dim < d_; dim++)
      {
        int compteur(0);
        //for(int k(0);k<resDonneesPartiel[i][dim].size();k++)
        for (vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
        {
          distRangPartiel[i][dim].push_back(distanceKendall(data_[dim][*it].x, resDonneesPartiel[i][dim][compteur]));
          compteur++;
        }
      }
    }
  }

  //changement de format
  vector<int> compteurElemPartiel(d_, 0);
  output_.distPartialRank = vector<vector<vector<int>>>(resDonneesPartiel.size());
  Rank rangTemp(d_);

  for (int iter(0); iter < (int)distRangPartiel.size(); iter++)
  {
    for (int dim(0); dim < d_; dim++)
      compteurElemPartiel[dim] = 0;

    for (int ind(0); ind < n_; ind++)
    {
      for (int dim(0); dim < d_; dim++)
      {
        if (data_[dim][ind].isPartial)
        {
          rangTemp[dim] = distRangPartiel[iter][dim][compteurElemPartiel[dim]];
          compteurElemPartiel[dim]++;
        }
        else
          rangTemp[dim] = 0;
      }
      output_.distPartialRank[iter].push_back(rangTemp);
    }
  }
}

void RankCluster::storeParameters(int const iterNumber, vector<vector<double>> &resProp,
                                  vector<vector<vector<double>>> &resP, vector<vector<vector<Rank>>> &resMu,
                                  vector<vector<int>> &resZ, vector<vector<vector<Rank>>> &resPartialData)
{
  // for identifiability, we must have p >= 0.5. If not, we invert the rank.
  for (int l(0); l < d_; l++)
  {
    for (int k(0); k < g_; k++)
    {
      if (p_[l][k] < 0.5)
      {
        p_[l][k] = 1 - p_[l][k];
        invertRank(mu_[l][k]);
      }
    }
  }

  // the first cluster must be the cluster with the fewer index of mu
  vector<int> indRank(g_);
  for (int k(0); k < g_; k++)
    indRank[k] = rank2index(mu_[0][k], tab_factorial(m_[0]));

  tri_insertionMulti(indRank, proportion_, p_, mu_, z_, g_, d_, n_);

  resP[iterNumber - parameter_.burnAlgo] = p_;
  resProp[iterNumber - parameter_.burnAlgo] = proportion_;
  resMu[iterNumber - parameter_.burnAlgo] = mu_;
  resZ[iterNumber - parameter_.burnAlgo] = z_;

  for (int dim(0); dim < d_; dim++)
  {
    int compteur(0);
    for (vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
    {
      resPartialData[iterNumber - parameter_.burnAlgo][dim][compteur] = data_[dim][*it].x;
      compteur++;
    }
  }

}

void RankCluster::run()
{
  convergence_ = false;
  int nbTry(0);
  while (!convergence_ && nbTry < parameter_.maxTry)
  {
    try
    {
      // double t0, t1, t2, t3, tM(0), tSE(0);

      // if (parameter_.verbose)
      // {
      //     cout << "##########################################################" << endl;
      //     cout << "#  SEM-Gibbs Algorithm for multivariate partial ranking  #" << endl;
      //     cout << "##########################################################" << endl;
      // }
      // t0 = clock();
      initialization();
      // t1 = clock();

      // if (parameter_.verbose)
      //     cout << "Initialization: " << (double)(t1 - t0) / CLOCKS_PER_SEC << "s." << endl;

      // objects for storing the estimated parameters at each iteration
      vector<vector<vector<double>>> resP(parameter_.maxIt - parameter_.burnAlgo,
                                          vector<vector<double>>(d_, vector<double>(g_)));
      vector<vector<double>> resProp(parameter_.maxIt - parameter_.burnAlgo, (vector<double>(g_)));
      vector<vector<int>> resZ(parameter_.maxIt - parameter_.burnAlgo, vector<int>(n_));
      vector<vector<vector<Rank>>> resMu(parameter_.maxIt - parameter_.burnAlgo, mu_);
      vector<vector<vector<Rank>>> resDonneesPartiel(parameter_.maxIt - parameter_.burnAlgo, output_.initialPartialRank);

      for (int iter(0); iter < parameter_.maxIt; iter++)
      {
        // if(parameter_.verbose)
        //    cout << "*";

        // t2 = clock();
        SEstep();
        // t3 = clock();
        // tSE += t3-t2;

        // t2 = clock();
        Mstep();
        // t3 = clock();
        // tM += t3-t2;

        // we store the estimated parameters
        if (iter >= parameter_.burnAlgo)
        {
          storeParameters(iter, resProp, resP, resMu, resZ, resDonneesPartiel);
        }
      }

      // t2 = clock();
      // if (parameter_.verbose)
      // {
      //     cout << endl << endl << "Log-likelihood estimation" << endl;

      //     cout << "Computing time for SE step: " << (double)tSE / CLOCKS_PER_SEC << "s ( "
      //          << (double)tSE / CLOCKS_PER_SEC / parameter_.maxIt << "s per step)." << endl;
      //     cout << "Computing time for M step: " << (double)tM / CLOCKS_PER_SEC << "s ( "
      //          << (double)tM / CLOCKS_PER_SEC / parameter_.maxIt << "s per step )." << endl;
      // }

      // compute log-likelihood and choice of the best parameters
      selectBestParameters(resMu, resP, resProp);
      //t3=clock();

      // compute the partition associated to the best result
      computePartition();

      // compute distance between estimated parameters at each iteration
      computeDistance(resProp, resP, resMu, resZ, resDonneesPartiel);

      // compute criterion
      output_.bic = BIC(output_.L, n_, 2 * g_ * d_ + g_ - 1);

      output_.icl = output_.bic;

      output_.entropy = ArrayXd(n_);
      for (int i = 0; i < n_; i++)
      {
        output_.entropy(i) = 0;
        for (int j = 0; j < g_; j++)
        {
          if (output_.tik(i, j) != 0)
            output_.entropy(i) -= output_.tik(i, j) * std::log(output_.tik(i, j));
        }
        output_.icl += 2 * output_.entropy(i);
      }

      // if (parameter_.verbose)
      // {
      //     cout << "Total computing time : " << (double)(t3 - t0) / CLOCKS_PER_SEC << "s" << endl;
      // }

      // if p < 0.5, we invert the associate mu and put p=1-p
      for (int j = 0; j < d_; j++)
        for (int k = 0; k < g_; k++)
        {
          if (p_[j][k] < 0.5)
          {
            p_[j][k] = 1 - p_[j][k];
            invertRank(mu_[j][k]);
          }
        }

        convergence_ = true;

    } // end try
    catch (string const &chaine)
    {
      convergence_ = false;
    }
    nbTry++;
  }
}

// ! not use in run but for interface with R
//return probability to belong to the final cluster
ArrayXd RankCluster::probability() const
{
  ArrayXd probability(n_);
  for (int ind(0); ind < n_; ind++)
    probability(ind) = output_.probabilities(ind, z_[ind]);

  return probability;
}

// ! not use in run but for interface with R
//LL gibbs
void RankCluster::estimateCriterion(double &L, double &bic, double &icl)
{

  /*initialisation partial rank and order of presentation*/
  //partial data and order of presentation initialization
  for (int dim(0); dim < d_; dim++)
  {
    Rank rankTemp(m_[dim]);
    initializeRank(rankTemp);
    for (int ind(0); ind < n_; ind++)
    {
      //initialization of y
      Rshuffle(rankTemp.begin(), rankTemp.end());
      data_[dim][ind].y = rankTemp;

      if (data_[dim][ind].isPartial)
      {
        for (int i = 0; i < (int)data_[dim][ind].missingData.size(); i++)
        {
          //initialization of Partial Rank
          Rank rankTemp2(data_[dim][ind].missingIndex[i]);
          Rshuffle(rankTemp2.begin(), rankTemp2.end());

          for (int ii = 0; ii < (int)data_[dim][ind].missingData[i].size(); ii++)
            data_[dim][ind].x[rankTemp2[ii]] = data_[dim][ind].missingData[i][ii];
        }
      }
    }
  }

  /*log likelihood computation*/
  ArrayXXd tik(n_, g_);
  long double p1(0), p2(0), p1x(0), p2x(0), alea(0), li(0);
  double div((double)(parameter_.nGibbsL - parameter_.burnL));
  vector<int> compteur(d_, 0);
  Rank x1, x2;

  //objet pour stockage de calcul pour éviter répétition
  ArrayXXd proba1(d_, g_), proba2(d_, g_);
  ArrayXd proba1X(g_), proba2X(g_);

  ArrayXd prop(g_);
  for (int i(0); i < g_; i++)
    prop(i) = proportion_[i];
  ArrayXXd propb(1, g_);
  propb = prop.transpose();

  //génération rang 1 2 3 ..m par dimension
  vector<Rank> yTemp(d_), x(d_);
  for (int j(0); j < d_; j++)
  {
    yTemp[j].resize(m_[j]);
    initializeRank(yTemp[j]);
  }

  vector<double> logL(parameter_.nGibbsL - parameter_.burnL, 0);

  //simulation de y multi dimensionnel
  for (int ind(0); ind < n_; ind++)
  {

    vector<Rank> y(d_), y2(d_), y1(d_);
    li = 0;
    //algorithme de Gibbs pour simuler yi

    //initialisation de y et p pour Gibbs
    y = yTemp;
    for (int j(0); j < d_; j++)
    {
      Rshuffle(y[j].begin(), y[j].end()); //permutation de 1 2 3 ..m
      x[j] = data_[j][ind].x;
    }

    y1 = y;

    for (int k(0); k < g_; k++)
    {
      tik(ind, k) = 0;
      for (int j(0); j < d_; j++)
        proba1(j, k) = probaCond(x[j], y1[j], mu_[j][k], p_[j][k]);
    }

    p1 = (long double)(propb * proba1.colwise().prod()).sum();
    proba2 = proba1;

    for (int iter(0); iter < parameter_.nGibbsL; iter++)
    {

      /*simulation des y*/
      for (int J(0); J < d_; J++)
      {
        for (int K(0); K < m_[J] - 1; K++)
        {
          //"état" à tester (inversion de 2 éléments adjacents)
          y2 = y;
          y2[J][K] = y[J][K + 1];
          y2[J][K + 1] = y[J][K];

          //tester un stockage des proba calculées pour éviter répétition de calculs dans la boucle
          for (int k(0); k < g_; k++)
            proba2(J, k) = probaCond(x[J], y2[J], mu_[J][k], p_[J][k]);

          p2 = (long double)(propb * proba2.colwise().prod()).sum();

          alea = (long double)runif(0., p1 + p2); //unif(0,p1+p2)

          if (alea < p2) //accept changement
          {
            y[J] = y2[J];
            p1 = p2;
            proba1.row(J) = proba2.row(J);
            y1[J] = y[J];
          }
          else //do not change y
          {
            y[J] = y1[J]; //rajout J
            proba2.row(J) = proba1.row(J);
          }
        }
      }
      //y_i is updated

      /*simulation of partial rank with a gibbs sampler*/
      for (int J = 0; J < d_; J++)
      {
        if (data_[J][ind].isPartial) //simulation of xi if it is a partial rank
        {
          x1 = x[J];
          proba1X = proba1.row(J);
          p1x = (proba1X * prop).sum();

          for (int kk = 0; kk < (int)(data_[J][ind].missingIndex).size() - 1; kk++) //Gibbs sur les x
          {
            for (int k = 0; k < (int)(data_[J][ind].missingIndex[kk]).size() - 1; k++)
            {
              //new x to test
              x2 = x[J];
              x2[data_[J][ind].missingIndex[kk][k]] = x[J][data_[J][ind].missingIndex[kk][k + 1]];
              x2[data_[J][ind].missingIndex[kk][k + 1]] = x[J][data_[J][ind].missingIndex[kk][k]];

              for (int l = 0; l < g_; l++)
                proba2X(l) = probaCond(x2, y[J], mu_[J][l], p_[J][l]);

              p2x = (proba2X * prop).sum();

              alea = (double)runif(0., p1x + p2x);

              if (alea < p2) //we accept the changement
              {
                x[J] = x2;
                p1x = p2x;
                proba1X = proba2X;
                x1 = x[J];
              }
              else
                x[J] = x1;
            }
          }
          proba1.row(J) = proba1X;
        }
      }

      if (iter >= parameter_.burnL)
      {
        ArrayXd calculInter(g_);
        for (int cl = 0; cl < g_; cl++)
        {
          calculInter(cl) = 1;
          for (int dim = 0; dim < d_; dim++)
            calculInter(cl) *= proba1(dim, cl);
          calculInter(cl) *= propb(cl);
        }

        double den = calculInter.sum();
        tik.row(ind) += (calculInter / den);

        li += (long double)1 / den;
      }

    } //end gibbs sampling for sample ind

    //L -= log(li*div);
    L -= log(li);

    tik.row(ind) /= div;

  } //end loop on sample

  L += (double)n_ * log(div);

  output_.L = L;

  output_.bic = BIC(output_.L, n_, 2 * g_ * d_ + g_ - 1);

  output_.icl = output_.bic;

  ArrayXd entropy(n_);
  for (int i = 0; i < n_; i++)
  {
    entropy(i) = 0;
    for (int j = 0; j < g_; j++)
    {
      if (tik(i, j) != 0)
        entropy(i) -= tik(i, j) * std::log(tik(i, j));
    }
    output_.icl += 2 * entropy(i);
  }

  bic = output_.bic;
  icl = output_.icl;
}
