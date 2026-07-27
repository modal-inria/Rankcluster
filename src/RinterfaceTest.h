#ifndef RUNTEST_H_
#define RUNTEST_H_

#include <RcppEigen.h>
#include <vector>
#include <set>
#include "test.h"

// Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
std::vector<std::vector<double>> convertToVVd(Rcpp::NumericMatrix data);

// Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
std::vector<std::vector<int>> convertToVVi(Rcpp::NumericMatrix data);

// Convert Rcpp::NumericMatrix to vector of ranks
std::vector<RankStruct> downUniVariateRank(Rcpp::NumericMatrix XR);

#endif /* RUNTEST_H_ */
