/*
 * File containing functions associated with the ISR distribution
 */
#include <map>

#include "ISRfunctions.h"
#include "functions.h"

using namespace std;

/*
 * Computation of A(x,y) (number of comparisons) and G(x,y,mu) (number of good comparisons)
 *
 * These values are useful to compute the probability according to an ISR.
 * See appendix A and B: https://hal.archives-ouvertes.fr/hal-00441209v3/document
 */
vector<int> comparaison(Rank const &x, Rank const &y, Rank const &mu)
{
  int const m(mu.size());
  int gplus(0), gmoins(0), gjmoinsb(0), gjplusb(0), index(0);
  vector<int> ajmoins, ajplus, ajplusb, ajmoinsb, ajplusbIndex;
  ajplusb.reserve(m);           //le Aj+ en cours
  ajmoinsb.reserve(m);          //le Aj- en cours
  ajplusbIndex.reserve(m);      //les index du Aj+ en cours
  ajplus.reserve(m * (m - 1));  //l'union de tt les Aj+
  ajmoins.reserve(m * (m - 1)); //l'union de tt les Aj-
  
  for (int j(1); j < m; j++)
  {
    gjmoinsb = 0;
    gjplusb = 0;
    for (int i(0); i < j; i++)
    {
      //calcul Aj-
      if (positionRank(x, y[i]) < positionRank(x, y[j]))
      {
        ajmoins.push_back(i);
        ajmoinsb.push_back(i);
      }
      else //calcul Aj+//if (positionRank(x,y[i])>positionRank(x,y[j]))
      {
        ajplusb.push_back(positionRank(x, y[i]));
        ajplusbIndex.push_back(i);
      }
    }
    
    if (ajplusb.size() > 0) //si le Aj+ en cours est non vide, on rajoute l'index du min à Aj+
    {
      index = min_element(ajplusb.begin(), ajplusb.end()) - ajplusb.begin();
      ajplus.push_back(ajplusbIndex[index]);
      
      //calcul de G+
      if (positionRank(mu, y[j]) < positionRank(mu, y[ajplus[ajplus.size() - 1]]))
      {
        gplus++;
        gjplusb++;
      }
      ajplusb.erase(ajplusb.begin(), ajplusb.end());
      ajplusbIndex.erase(ajplusbIndex.begin(), ajplusbIndex.end());
    }
    if (ajmoinsb.size() > 0) //si le Aj- en cours est non vide on calcule G-
    {
      //calcul de G-
      for (unsigned int i(0); i < ajmoinsb.size(); i++)
      {
        if (positionRank(mu, y[ajmoinsb[i]]) < positionRank(mu, y[j]))
        {
          gmoins++;
          gjmoinsb++;
        }
      }
      ajmoinsb.erase(ajmoinsb.begin(), ajmoinsb.end());
    }
  }
  
  vector<int> comparaison(2, 0);
  comparaison[0] = ajmoins.size() + ajplus.size();
  comparaison[1] = gmoins + gplus;
  
  return comparaison;
}

/** Compute log p(x|y;mu,pi) according to an ISR model */
double lnProbaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p)
{
  vector<int> comp(2);
  comp = comparaison(x, y, mu);
  
  // manage special case, otherwise it returns NaN
  if (p == 1)
  {
    if (comp[0] == comp[1])
      return 0.;
  }
  
  if (p == 0)
  {
    if (comp[1] == 0)
      return 0.;
  }
  
  return comp[1] * log(p) + (comp[0] - comp[1]) * log(1. - p);
}

double probaCond(Rank const &x, Rank const &y, Rank const &mu, double const &p)
{
  return exp(lnProbaCond(x, y, mu, p));
}


// Simulate n samples from an ISR(mu, p)
vector<Rank> simulISR(int const &n, int const &m, Rank const &mu, double const &p)
{
  vector<Rank> simul(n, Rank(m, 0));
  Rank s(m), rgTemp(m);
  int l;
  double correct;
  bool compar, avance;
  
  initializeRank(rgTemp);
  
  for (int i(0); i < n; i++)
  {
    //simulation d'un rang aléatoire: permutation du vecteur 1 2..m
    s = rgTemp;
    Rshuffle(s.begin(), s.end());
    simul[i][0] = s[0];
    for (int j(1); j < m; j++)
    {
      l = 0;
      avance = true;
      while (avance && l < j)
      {
        correct = (double)runif(0., 1.);
        compar = (positionRank(mu, s[j]) < positionRank(mu, simul[i][l]));
        if ((compar && correct < p) || (!compar && correct > p))
        {
          for (int k(j - 1); k >= l; k--)
            simul[i][k + 1] = simul[i][k];
          
          simul[i][l] = s[j];
          avance = false;
        }
        else
          l++;
      }
      if (l == j)
        simul[i][l] = s[j];
    }
  }
  return simul;
}


// Simulate n samples from a mixture of ISR(mu, p)
void simulMixtureISR(vector<Rank> &simul, vector<Rank> const &mu,
                     vector<double> const &p, vector<double> const &prop)
{
  int classe = 0;
  int n = simul.size();
  int m = mu[0].size();
  
  Eigen::ArrayXd probaEig(prop.size());
  for (size_t k = 0; k < prop.size(); k++)
    probaEig(k) = prop[k];
  
  for (int i(0); i < n; i++)
  {
    // class number sample according to proportion
    classe = sampleMultinomial(probaEig);
    
    // sample a rank according to the parameters of the sampled class
    simul[i] = simulISR(1, m, mu[classe], p[classe])[0];
  }
}

double proba(vector<Rank> const &x, vector<Rank> const &mu, vector<double> const &pi)
{
  int d = pi.size();
  vector<double> probaDim(d, 0);
  double finalProba;
  vector<int> tabFact, listeY;
  Rank y;
  vector<int> m(mu.size());
  map<int, vector<int>> diffDim;
  
  for (int dim = 0; dim < d; dim++)
  {
    m[dim] = mu[dim].size();
    diffDim[m[dim]].push_back(dim);
  }
  
  // we iterate over the different values of m (size of ranks) to limit the computation of y
  for (map<int, vector<int>>::iterator it = diffDim.begin(); it != diffDim.end(); it++)
  {
    int m = it->first;
    tabFact = tab_factorial(m);
    listeY = listIndexOrderOfPresentation(m, tabFact);
    double mult = 2. / (double)tabFact[m - 1];
    
    for (int i = 0; i < tabFact[m - 1] / 2; i++)
    {
      y = index2rankNoCheck(listeY[i], m, tabFact);
      
      // for each dimension that have a rank of size m
      for (int j = 0; j < (int)(it->second).size(); j++)
        probaDim[(it->second)[j]] += probaCond(x[(it->second)[j]], y, mu[(it->second)[j]], pi[(it->second)[j]]);
    }
    
    for (int j = 0; j < (int)(it->second).size(); j++)
      probaDim[(it->second)[j]] *= mult;
  }
  
  finalProba = probaDim[0];
  for (int dim = 1; dim < d; dim++)
    finalProba *= probaDim[dim];
  
  return finalProba;
}
