/*
 * rank_utils.cpp
 *
 * Utilities for rank manipulation.
 */

#include <algorithm>  // std::find, std::remove, std::swap, std::min_element
#include <numeric>    // std::iota

#include "rank_utils.h"
#include "sampling_utils.h"   // Rshuffle (used by randomRank)

int positionRank(Rank const &x, int const &i)
{
    int j = 0;
    while (x[j] != i)
        j++;
    return j;
}

void initializeRank(Rank &rank)
{
    std::iota(rank.begin(), rank.end(), 1);
}

void randomRank(Rank &rank)
{
    initializeRank(rank);
    Rshuffle(rank.begin(), rank.end());
}

int factorial(int const &nombre)
{
    if (nombre <= 1)
        return 1;
    return nombre * factorial(nombre - 1);
}

std::vector<int> tab_factorial(int const &m)
{
    std::vector<int> tab(m);
    tab[0] = 1;
    for (int i = 1; i < m; i++)
        tab[i] = (i + 1) * tab[i - 1];
    return tab;
}

int rank2index(Rank const &rank, std::vector<int> const &tabFactorial)
{
    int const m(static_cast<int>(rank.size()));
    int index = (rank[0] - 1) * tabFactorial[m - 2];

    std::vector<int> liste(m);
    initializeRank(liste);
    liste.erase(std::remove(liste.begin(), liste.end(), rank[0]), liste.end());

    for (int j = 1; j < m - 1; j++)
    {
        auto foundIt = std::find(liste.begin(), liste.end(), rank[j]);
        index += static_cast<int>(foundIt - liste.begin()) * tabFactorial[m - j - 2];
        liste.erase(std::remove(liste.begin(), liste.end(), rank[j]), liste.end());
    }

    return index + 1;
}

std::vector<int> rank2index(std::vector<Rank> const &rankList,
                              std::vector<int> const &tabFact)
{
    int n(static_cast<int>(rankList.size()));
    std::vector<int> listeIndex(n);
    for (int i = 0; i < n; i++)
        listeIndex[i] = rank2index(rankList[i], tabFact);
    return listeIndex;
}

Rank index2rankNoCheck(int index, int const &m,
                        std::vector<int> const &tabFactorial)
{
    Rank r(m, 0);
    std::vector<int> liste(m);

    index--;
    r[0] = index / tabFactorial[m - 2] + 1;

    initializeRank(liste);
    liste.erase(std::remove(liste.begin(), liste.end(), r[0]), liste.end());

    for (int j = 1; j < m - 1; j++)
    {
        int temp = index;
        for (int k = 1; k < j + 1; k++)
            temp %= tabFactorial[m - k - 1];

        int temp2 = temp / tabFactorial[m - j - 2];
        r[j] = liste[temp2];
        liste.erase(std::remove(liste.begin(), liste.end(), r[j]), liste.end());
    }
    r[m - 1] = liste[0];

    return r;
}

Rank index2rank(int index, int const &m, std::vector<int> const &tabFactorial)
{
    if (index > factorial(m))
        return Rank(m, 0);
    return index2rankNoCheck(index, m, tabFactorial);
}

Rank index2rank(int index, int const &m)
{
    std::vector<int> tabFactorial = tab_factorial(m);
    return index2rank(index, m, tabFactorial);
}

std::vector<int> listIndexOrderOfPresentation(int const &m,
                                               std::vector<int> const &tabFactorial)
{
    std::vector<int> liste(tabFactorial[m - 1] / 2);
    int startIdx = 0, endIdx = 0, ind = 0;

    for (int i = 1; i <= m - 1; i++)
    {
        startIdx = (i - 1) * (tabFactorial[m - 2] + tabFactorial[m - 3]) + 1;
        endIdx   = i * tabFactorial[m - 2];
        int nb   = endIdx - startIdx + 1;
        for (int j = 0; j < nb; j++)
            liste[ind + j] = startIdx + j;
        ind += nb;
    }

    return liste;
}

void invertRank(Rank &rank)
{
    int n = static_cast<int>(rank.size());
    for (int j = 0; j < n / 2; j++)
        std::swap(rank[j], rank[n - j - 1]);
}

Rank ordering2ranking(Rank const &x)
{
    Rank y(x);
    for (int i = 0; i < static_cast<int>(x.size()); i++)
        y[i] = positionRank(x, i + 1) + 1;
    return y;
}

int distanceKendall(Rank const &x, Rank const &y)
{
    int const m(static_cast<int>(x.size()));
    Rank xr = ordering2ranking(x);
    Rank yr = ordering2ranking(y);
    int distK = 0;

    for (int i = 0; i < m - 1; i++)
        for (int j = i + 1; j < m; j++)
            if ((xr[i] - xr[j]) * (yr[i] - yr[j]) < 0)
                distK++;

    return distK;
}
