/*
 * functions.cpp
 *
 *  Created on: 1 mars 2013
 *      Author: Quentin Grimonprez
 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <Rmath.h>

#include "functions.h"

using namespace std;

// generate an integer between 0 and n -1
int randWrapper(const int n)
{
    return floor(unif_rand() * n);
}

// generate a rank containing 1, 2, .., m
void initializeRank(vector<int> &rank)
{
    for (size_t j = 0; j < rank.size(); j++)
        rank[j] = j + 1;
}

// generate a random Rank
void randomRank(vector<int> &rank)
{
    initializeRank(rank);
    Rshuffle(rank.begin(), rank.end());
}

// return the position of object i in x (ordering representation). Equivalent rank.gamma
// TODO manage the case where is is not in x ?
int positionRank(vector<int> const &x, int const &i)
{
    int j(0);
    while (x[j] != i)
        j++;

    return j; //entre 0 et m-1
}


// factorial function
int factorial(int const &m)
{
    int temp;

    if (m <= 1)
        return 1;

    temp = m * factorial(m - 1);
    return temp;
}

// compute all factorial values from 1! to to m!
vector<int> tab_factorial(int const &m)
{
    vector<int> tab(m);
    tab[0] = 1;
    for (int i(1); i < m; i++)
        tab[i] = (i + 1) * tab[i - 1];

    return tab;
}

//conversion d'un rang en index
int rank2index(vector<int> const &rang, vector<int> const &tabFactorial)
{
    int const m(rang.size());
    int index(0);
    index = (rang[0] - 1) * tabFactorial[m - 2];
    vector<int>::iterator it;
    vector<int> liste(m);

    initializeRank(liste);

    liste.erase(remove(liste.begin(), liste.end(), rang[0]), liste.end());

    for (int j(1); j < m - 1; j++)
    {
        it = search_n(liste.begin(), liste.end(), 1, rang[j]);
        index += (int(it - liste.begin()) * tabFactorial[m - j - 2]);
        liste.erase(remove(liste.begin(), liste.end(), rang[j]), liste.end());
    }

    return index + 1;
}

vector<int> rank2index(vector<vector<int>> const &listeRang, vector<int> const &tabFact)
{
    const unsigned int n(listeRang.size());
    vector<int> listeIndex(n);
    for (unsigned int i(0); i < n; i++)
        listeIndex[i] = rank2index(listeRang[i], tabFact);

    return listeIndex;
}

//conversion d'un index en rang
vector<int> index2rank(int index, int const &m, vector<int> const &tabFactorial)
{
    if (index > factorial(m))
    {
        vector<int> temp(m, 0);
        return temp;
        //cout<<"ERREUR "<<index<<"<"<<m<<"!"<<endl;
    }
    else
    {
        vector<int> liste(m), r(m, 0);
        int temp(0), temp2(0);

        index--;
        r[0] = index / tabFactorial[m - 2] + 1;

        initializeRank(liste);

        //on supprime l'élément égale à r[0]
        liste.erase(remove(liste.begin(), liste.end(), r[0]), liste.end());

        for (int j(1); j < m - 1; j++)
        {
            temp = index;
            for (int k(1); k < j + 1; k++)
                temp %= tabFactorial[m - k - 1];

            temp2 = temp / tabFactorial[m - j - 2];
            r[j] = liste[temp2];

            //replace (liste.begin(), liste.end(), r[j], 0);
            liste.erase(remove(liste.begin(), liste.end(), r[j]), liste.end());
        }
        r[m - 1] = liste[0];

        return r;
    }
}

vector<int> index2rank(int index, int const &m)
{
    vector<int> r(m, 0);

    if (index > factorial(m))
    {
        vector<int> temp(m, 0);
        return temp;
        //cout<<"ERREUR "<<index<<"<"<<m<<"!"<<endl;
    }
    else
    {
        vector<int> liste(m);
        int temp(0), temp2(0);

        index--;
        r[0] = index / factorial(m - 1) + 1;

        initializeRank(liste);

        //on supprime l'élément égale à r[0]
        liste.erase(remove(liste.begin(), liste.end(), r[0]), liste.end());

        for (int j(1); j < m - 1; j++)
        {
            temp = index;
            for (int k(1); k < j + 1; k++)
                temp %= factorial(m - k);

            temp2 = temp / factorial(m - j - 1);
            r[j] = liste[temp2];

            //replace (liste.begin(), liste.end(), r[j], 0);
            liste.erase(remove(liste.begin(), liste.end(), r[j]), liste.end());
        }
        r[m - 1] = liste[0];

        return r;
    }
}

vector<int> index2rankb(int index, int const &m, vector<int> const &tabFactorial)
{
    vector<int> liste(m), r(m, 0);
    int temp(0), temp2(0);

    index--;
    r[0] = index / tabFactorial[m - 2] + 1;

    initializeRank(liste);

    //on supprime l'élément égale à r[0]
    liste.erase(remove(liste.begin(), liste.end(), r[0]), liste.end());

    for (int j(1); j < m - 1; j++)
    {
        temp = index;
        for (int k(1); k < j + 1; k++)
            temp %= tabFactorial[m - k - 1];

        temp2 = temp / tabFactorial[m - j - 2];
        r[j] = liste[temp2];

        //replace (liste.begin(), liste.end(), r[j], 0);
        liste.erase(remove(liste.begin(), liste.end(), r[j]), liste.end());
    }
    r[m - 1] = liste[0];

    return r;
}

vector<int> listeSigma(int const &m, vector<int> const &tabFactorial)
{
    vector<int> liste(tabFactorial[m - 1] / 2);
    int const1(0), const2(0), ind(0);

    for (int i(1); i <= m - 1; i++)
    {
        const1 = (i - 1) * (tabFactorial[m - 2] + tabFactorial[m - 3]) + 1;
        const2 = i * tabFactorial[m - 2];
        int nb(const2 - const1 + 1);
        for (int j(0); j < nb; j++)
            liste[ind + j] = const1 + j;
        ind += nb;
    }

    return liste;
}

//inversion d'un rang
//ex: transforme 1 2 3 4 en 4 3 2 1
void inverseRang(vector<int> &rang)
{
    for (unsigned int j(0); j < rang.size() / 2; j++)
        swap(rang[j], rang[rang.size() - j - 1]);
}

void inverseRang(vector<int> &rang, int const &m)
{
    for (int j(0); j < m / 2; j++)
        swap(rang[j], rang[m - j - 1]);
}

//calcul BIC
double BIC(double loglikelihood, int nbDonnees, int nbParam)
{
    double bic(0);
    bic = -2 * loglikelihood + nbParam * log(nbDonnees);
    return bic;
}

double indiceRand(vector<int> const &z1, vector<int> const &z2)
{
    const int N(z1.size());
    double a(0), b(0), c(0), d(0);
    for (int i(0); i < N; i++)
    {
        for (int j(0); j < N; j++)
        {
            if ((z1[i] == z1[j]) && (z2[i] == z2[j]))
                a++;
            else
            {
                if ((z1[i] != z1[j]) && (z2[i] != z2[j]))
                    b++;
                else
                {
                    if ((z1[i] == z1[j]) && (z2[i] != z2[j]))
                        c++;
                    else
                        d++;
                }
            }
        }
    }

    return (a + b) / (a + b + c + d);
}

/*vector<int> order2rank(vector<int> &x,int const& m)
{
vector<int> y(x);
for(int i(0);i<m;i++)
x[i]=positionRank(y,i+1)+1;

return x;
}*/

vector<int> order2rank(vector<int> const &x, int const &m)
{
    vector<int> y(x);
    for (int i(0); i < m; i++)
        y[i] = positionRank(x, i + 1) + 1;

    return y;
}

int distanceKendall(vector<int> const &x, vector<int> const &y)
{
    const int m(x.size());
    vector<int> xr(m), yr(m);
    xr = order2rank(x, m);
    yr = order2rank(y, m);
    int distK(0);

    for (int i(0); i < m - 1; i++)
        for (int j(i + 1); j < m; j++)
            if ((xr[i] - xr[j]) * (yr[i] - yr[j]) < 0)
                distK++;

    return distK;
}

// mu: index des elements de la 1ere dim
//listeMu: listMu[dim][composante][elem]
void tri_insertionMulti(vector<int> &mu, vector<double> &prop, vector<vector<double>> &p, vector<vector<vector<int>>> &listeMu, vector<int> &z, int const &g, int const &d, int const &n)
{
    int i, j, elem;
    double elemprop;
    vector<double> elemp(d);
    vector<vector<int>> elemmu(d);
    vector<int> order(g);
    for (int l = 0; l < g; l++)
        order[l] = l;
    int elemorder;

    //sort algorithm
    //we sort the cluster of the first dim and make changement for all dim
    for (i = 1; i < g; ++i)
    {
        elem = mu[i];
        for (int l(0); l < d; l++)
            elemp[l] = p[l][i];
        elemprop = prop[i];
        for (int l(0); l < d; l++)
            elemmu[l] = listeMu[l][i];
        elemorder = order[i];

        for (j = i; j > 0 && mu[j - 1] > elem; j--)
        {
            order[j] = order[j - 1];
            mu[j] = mu[j - 1];
            prop[j] = prop[j - 1];
            for (int l(0); l < d; l++)
                p[l][j] = p[l][j - 1];
            for (int l(0); l < d; l++)
                listeMu[l][j] = listeMu[l][j - 1];
        }
        order[j] = elemorder;
        mu[j] = elem;
        for (int l(0); l < d; l++)
            p[l][j] = elemp[l];
        prop[j] = elemprop;

        for (int l(0); l < d; l++)
            listeMu[l][j] = elemmu[l];
    }

    //re order the z
    for (int l = 0; l < n; l++)
    {
        for (int k = 0; k < g; k++)
        {
            if (z[l] == order[k])
            {
                z[l] = k;
                break;
            }
        }
    }
}


//frequence d'un jeu de données multi varié
typedef struct TableauRang TableauRang;
pair<vector<vector<vector<int>>>, vector<int>> freqMulti(vector<vector<vector<int>>> const &listeRang)
{
    //rassemblement des mu identiques, proba->moyenne des proba,p->moyenne des p
    struct TableauRang;
    struct TableauRang
    {
        int compteur;
        std::vector<std::vector<int>> rang;
        TableauRang *suivant;
    };

    bool continuer(true), egaliteRang;
    const int d(listeRang.size()), N(listeRang[0].size());
    vector<vector<int>> rang(d);
    vector<int> M(d);
    for (int i(0); i < d; i++)
        M[i] = listeRang[i][0].size();

    TableauRang *headRang = new TableauRang;
    TableauRang *currRang = headRang;
    TableauRang *next = 0;
    currRang->compteur = 1;
    currRang->suivant = 0;
    for (int i(0); i < d; i++)
        rang[i] = listeRang[i][0];
    currRang->rang = rang;

    for (int j(1); j < N; j++)
    {
        continuer = true;
        currRang = headRang;
        while (continuer)
        {
            egaliteRang = true;
            for (int J(0); J < d; J++)
            {
                for (int k(0); k < M[J]; k++)
                {
                    if (currRang->rang[J][k] != listeRang[J][j][k])
                    {
                        egaliteRang = false;
                        break;
                    }
                }
            }

            if (egaliteRang)
            {
                currRang->compteur++;
                continuer = false;
            }
            else
            {
                if (currRang->suivant == 0)
                {
                    continuer = false;
                    next = new TableauRang;
                    currRang->suivant = next;
                    currRang = next;
                    currRang->compteur = 1;
                    currRang->suivant = 0;
                    for (int i(0); i < d; i++)
                        rang[i] = listeRang[i][j];
                    currRang->rang = rang;
                }
                else
                    currRang = currRang->suivant;
            }
        }
    }

    //stockage des infos
    vector<vector<vector<int>>> donnees(d);
    vector<int> freq;

    currRang = headRang;
    while (currRang->suivant != 0)
    {
        for (int i(0); i < d; i++)
            donnees[i].push_back(currRang->rang[i]);

        freq.push_back(currRang->compteur);
        next = currRang->suivant;
        delete currRang;
        currRang = next;
    }

    for (int i(0); i < d; i++)
        donnees[i].push_back(currRang->rang[i]);

    freq.push_back(currRang->compteur);
    delete currRang;

    return make_pair(donnees, freq);
}


bool acceptChange(double const logP1, double const logP2)
{
    double logP1PlusP2;
    if (logP1 > logP2)
    {
        logP1PlusP2 = logP1 + log(exp(logP2 - logP1) + 1);
    }
    else
    {
        logP1PlusP2 = logP2 + log(exp(logP1 - logP2) + 1);
    }

    double proba1 = exp(logP1 - logP1PlusP2);
    double proba2 = exp(logP2 - logP1PlusP2);

    double alea = runif(0., proba1 + proba2);

    return (alea < proba2);
}
