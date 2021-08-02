#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <vector>

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

struct ListeMu
{
    double compteur; //number of same mu
    std::vector<std::vector<Rank>> rangComplet;
    std::vector<std::vector<double>> p;
    std::vector<double> prop;
    ListeMu *suivant;
};


#endif // TYPEDEF_H